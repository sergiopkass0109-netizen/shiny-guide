/*
 * parallel.js — multi-core exact prime counting: the verified Lucy_Hedgehog
 * recurrence over SharedArrayBuffer, K threads, no locks in the hot loops.
 *
 * Schedule.  Primes are processed exactly as in engine.c: small primes one
 * at a time, larger primes in blocks of up to 32 consecutive primes
 * p1 < … < pk ≤ 2·p1.  For a block, every update of large[i] with
 * i > I0 = max(r/p1, x/p1³) reads only small[q] with q < p1² — a region no
 * small-pass of the block touches — so all k passes over large[(I0, imax]]
 * commute and are applied in ONE cache-blocked sweep (segments of SEG
 * entries, each prime's runs into a per-thread difference array, one
 * prefix-sum pass per segment).  The sweep has no data hazards at all: the
 * threads simply split the index range.  The prefix i ≤ I0 is then done
 * prime by prime in the classical order; there, and for the small primes,
 * the pass writes large[i] while reading large[i·p] (a HIGHER index) and
 * writes small[v] while reading small[⌊v/p⌋] (a LOWER index), so geometric
 * bands — [p^k, p^{k+1}) ascending for `large`, (r/p^{k+1}, r/p^k] descending
 * for `small` — read only indices outside themselves that no band has
 * written yet in this pass, and each band can be split across threads
 * freely.  One barrier per band or sweep.  Primes with little work
 * (imax < 2^16) are done by the coordinator alone, exactly as in the
 * single-thread engine.  Results are cross-checked against the three
 * single-thread engines in the test suite.
 *
 * Runs in Node (worker_threads) and in browsers that are cross-origin
 * isolated (the app's service worker arranges that on GitHub Pages).  The
 * same source runs in every thread; thread 0 coordinates.  Kernels are the
 * shared-memory WebAssembly build of engine_par.c; every kernel has a
 * line-for-line JavaScript twin used where that module cannot load.
 */
(function (root, factory) {
  "use strict";
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.NthPrimeParallel = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  var IS_NODE = typeof process !== "undefined" && !!(process.versions && process.versions.node) &&
                typeof require === "function";
  var KMAX = 32;                 // primes per block (matches engine.c)
  var BP_LEN = 8 + 5 * KMAX;     // block-parameter table: header + 5 doubles per prime

  // ---- thread body: self-contained, serialised with toString() ----------
  function workerMain() {
    var GEN = 0, REM = 1, NEXT = 2;   // control words: generation, threads remaining, next chunk
    var T_INIT = 1, T_PREFIX = 2, T_SMALL = 3, T_SWEEP = 4, T_EXIT = 9;
    var PAR_MIN = 1 << 16;   // primes with less work than this: coordinator alone
    var SEQ_MIN = 1 << 12;   // band ranges smaller than this: not worth waking the threads
    var DYN_MIN = 1 << 10;   // dynamically shared ranges are worth a barrier from here on
    var KMAX = 32, BP_LEN = 8 + 5 * KMAX;
    var isNode = typeof self === "undefined";
    var port = isNode ? require("worker_threads").parentPort : null;
    function post(m) { if (isNode) port.postMessage(m); else self.postMessage(m); }
    function onMessage(f) { if (isNode) port.on("message", f); else self.onmessage = function (e) { f(e.data); }; }

    onMessage(function (m) {
      var K = m.K, id = m.id, x = m.x, r = m.r, SEG = m.seg || 4096;
      var ctrl, par, bp, small, large, D, kern = null;
      if (m.memory) {
        // one shared WebAssembly memory: typed-array views for JS, offsets for the kernels
        var buf = m.memory.buffer;
        ctrl = new Int32Array(buf, m.ctrlOff, 16);
        par = new Float64Array(buf, m.parOff, 8);
        bp = new Float64Array(buf, m.bpOff, BP_LEN);
        small = new Uint32Array(buf, m.smallOff, r + 2);
        large = new Float64Array(buf, m.largeOff, r + 2);
        D = new Float64Array(buf, m.dOff + id * 8 * (SEG + 2), SEG + 2);
        try {
          var bin = m.parB64, bytes;
          if (typeof atob === "function") { var str = atob(bin); bytes = new Uint8Array(str.length); for (var bi = 0; bi < str.length; bi++) bytes[bi] = str.charCodeAt(bi); }
          else bytes = new Uint8Array(Buffer.from(bin, "base64"));
          var inst = new WebAssembly.Instance(new WebAssembly.Module(bytes), { env: { memory: m.memory } });
          kern = inst.exports;
        } catch (e) { kern = null; }
      } else {
        ctrl = new Int32Array(m.ctrl);
        par = new Float64Array(m.par);
        bp = new Float64Array(m.bp);
        small = new Uint32Array(m.small);
        large = new Float64Array(m.large);
        D = new Float64Array(SEG + 2);
      }
      var SO = m.smallOff || 0, LO = m.largeOff || 0, BP = m.bpOff || 0, CT = m.ctrlOff || 0;
      var DO = (m.dOff || 0) + id * 8 * (SEG + 2), X = BigInt(x);

      // this thread's slice of the index range [lo, hi]
      function share(lo, hi) {
        var len = hi - lo + 1;
        if (len <= 0) return null;
        var per = Math.ceil(len / K);
        var a = lo + id * per;
        var b = Math.min(hi, a + per - 1);
        return a <= b ? [a, b] : null;
      }

      function isqrtJ(n) {
        var s = Math.floor(Math.sqrt(n));
        while ((s + 1) * (s + 1) <= n) s++;
        while (s * s > n) s--;
        return s;
      }

      // ---- JavaScript twins of the engine_par.c kernels ----
      // large[i] -= small[⌊xp/i⌋] - sp1 for i in [a, b]; g = isqrt(xp).  Head
      // (i < g): one division each; then runs of equal quotient q, walked
      // downward — the run of q is (⌊xp/(q+1)⌋, ⌊xp/q⌋].
      function tailJS(a, b, xp, sp1, g) {
        var i = a, h = Math.min(b + 1, Math.max(a, g));
        for (; i < h; i++) large[i] -= small[Math.floor(xp / i)] - sp1;
        if (i > b) return;
        var q = Math.floor(xp / i), qMin = Math.max(1, Math.floor(xp / b)), ePrev = i - 1;
        for (; q >= qMin; q--) {
          var e = Math.min(b, Math.floor(xp / q)), c = small[q] - sp1;
          for (var j = ePrev + 1; j <= e; j++) large[j] -= c;
          ePrev = e;
        }
      }
      // sweep twin: runs go to the difference array D (indexed i - base)
      function tailDJS(a, b, xp, sp1, g, base) {
        var i = a, h = Math.min(b + 1, Math.max(a, g));
        for (; i < h; i++) large[i] -= small[Math.floor(xp / i)] - sp1;
        if (i > b) return;
        var q = Math.floor(xp / i), qMin = Math.max(1, Math.floor(xp / b)), ePrev = i - 1;
        for (; q >= qMin; q--) {
          var e = Math.min(b, Math.floor(xp / q)), c = small[q] - sp1;
          D[ePrev + 1 - base] += c;
          D[e + 1 - base] -= c;
          ePrev = e;
        }
      }
      // prefix pass of block prime j over [a, b] (b ≤ min(I0, IMAX[j]))
      function prefixJS(j, a, b) {
        var k = bp[0], I0 = bp[1];
        var p = bp[8 + 5 * j], xp = bp[9 + 5 * j], sp1 = bp[10 + 5 * j], g = bp[12 + 5 * j];
        var isw = Math.min(b, Math.floor(r / p)), i1 = Math.min(isw, Math.floor(I0 / p));
        var i = a;
        for (; i <= i1; i++) large[i] -= large[i * p] - sp1;
        for (; i <= isw; i++) {
          // I0 < m ≤ r: add back the sweep's contributions of primes l ≥ j
          var mm = i * p, val = large[mm];
          for (var l = j; l < k && bp[11 + 5 * l] >= mm; l++) val += small[Math.floor(bp[9 + 5 * l] / mm)] - bp[10 + 5 * l];
          large[i] -= val - sp1;
        }
        tailJS(i, b, xp, sp1, g);
      }
      // cache-blocked sweep of [a, b] (a > I0, b ≤ IMAX[0]) for every prime of the block
      function sweepJS(a, b) {
        var k = bp[0], t, j;
        for (var L = a; L <= b;) {
          var R = Math.min(b, L + SEG - 1), n = R - L + 1;
          for (t = 0; t <= n; t++) D[t] = 0;
          for (j = 0; j < k && bp[11 + 5 * j] >= L; j++) tailDJS(L, Math.min(R, bp[11 + 5 * j]), bp[9 + 5 * j], bp[10 + 5 * j], bp[12 + 5 * j], L);
          var acc = 0;
          for (t = 0; t < n; t++) { acc += D[t]; large[L + t] -= acc; }
          L = R + 1;
        }
      }
      // dynamic sharing of [lo, hi]: chunk s is [lo + s·ch, lo + (s+1)·ch − 1]
      function sweepDynJS(lo, hi, ch) {
        for (;;) {
          var L = lo + Atomics.add(ctrl, NEXT, 1) * ch;
          if (L > hi) break;
          sweepJS(L, Math.min(hi, L + ch - 1));
        }
      }
      function prefixDynJS(j, lo, hi, ch) {
        for (;;) {
          var L = lo + Atomics.add(ctrl, NEXT, 1) * ch;
          if (L > hi) break;
          prefixJS(j, L, Math.min(hi, L + ch - 1));
        }
      }
      function blockTailJS() {
        var k = bp[0], I0 = bp[1];
        for (var j = 0; j < k; j++) {
          var imaxj = Math.min(I0, bp[11 + 5 * j]);
          if (imaxj >= 1) prefixJS(j, 1, imaxj);
          var p = bp[8 + 5 * j], p2 = p * p;
          if (p2 <= r) smallJS(p2, r, p, bp[10 + 5 * j]);
        }
      }
      function smallJS(a, b, p, sp1) {
        for (var v = b; v >= a;) {
          var q = Math.floor(v / p), w = Math.max(q * p, a), sub = small[q] - sp1;
          for (; v >= w; v--) small[v] -= sub;
        }
      }

      function task(type) {
        var c, v;
        if (type === T_INIT) {
          c = share(1, r);
          if (!c) return;
          if (kern) { kern.init_range(SO, LO, X, BigInt(c[0]), BigInt(c[1])); return; }
          for (v = c[0]; v <= c[1]; v++) { small[v] = v - 1; large[v] = Math.floor(x / v) - 1; }
        } else if (type === T_PREFIX) {
          if (par[4]) {                                   // hazard-free tail: dynamic chunks of par[4]
            if (kern) kern.prefix_dyn(SO, LO, BP, CT, par[1], BigInt(par[2]), BigInt(par[3]), BigInt(par[4]));
            else prefixDynJS(par[1], par[2], par[3], par[4]);
            return;
          }
          c = share(par[2], par[3]);                      // band: equal shares, uniform cost
          if (!c) return;
          if (kern) kern.prefix_range(SO, LO, BP, par[1], BigInt(c[0]), BigInt(c[1]));
          else prefixJS(par[1], c[0], c[1]);
        } else if (type === T_SWEEP) {
          if (kern) kern.sweep_dyn(SO, LO, BP, DO, CT, BigInt(par[1]), BigInt(par[2]), BigInt(par[3]));
          else sweepDynJS(par[1], par[2], par[3]);
        } else if (type === T_SMALL) {
          c = share(par[1], par[2]);
          if (!c) return;
          if (kern) kern.small_range(SO, BigInt(c[0]), BigInt(c[1]), BigInt(par[3]), BigInt(par[4]));
          else smallJS(c[0], c[1], par[3], par[4]);
        }
      }

      if (id !== 0) {
        var seen = 0;
        for (;;) {
          Atomics.wait(ctrl, GEN, seen);
          seen = Atomics.load(ctrl, GEN);
          var type = par[0];
          if (type === T_EXIT) break;
          task(type);
          Atomics.sub(ctrl, REM, 1);
          Atomics.notify(ctrl, REM);
        }
        return;
      }

      // ---- coordinator (thread 0) ----
      function dispatch(type, a1, a2, a3, a4) {
        par[0] = type; par[1] = a1; par[2] = a2; par[3] = a3; par[4] = a4;
        Atomics.store(ctrl, NEXT, 0);
        Atomics.store(ctrl, REM, K - 1);
        Atomics.add(ctrl, GEN, 1);
        Atomics.notify(ctrl, GEN);
        task(type);
        var rem;
        while ((rem = Atomics.load(ctrl, REM)) > 0) Atomics.wait(ctrl, REM, rem);
      }
      // chunk size for a dynamically shared range: a few chunks per thread, never above a segment
      function chunkFor(len) { return Math.max(256, Math.min(SEG, Math.ceil(len / (4 * K)))); }
      function prefixRange(j, a, b) {
        if (a > b) return;
        if (kern) kern.prefix_range(SO, LO, BP, j, BigInt(a), BigInt(b));
        else prefixJS(j, a, b);
      }

      // block parameters (shared with every thread and both kernel kinds)
      var P = new Array(KMAX);
      function setBlock(k, I0) {
        bp[0] = k; bp[1] = I0; bp[2] = r;
        for (var j = 0; j < k; j++) {
          var pj = P[j], xp = Math.floor(x / pj);
          bp[8 + 5 * j] = pj;
          bp[9 + 5 * j] = xp;
          bp[10 + 5 * j] = small[pj - 1];                       // π(p_j − 1): final, p_j − 1 < p1²
          bp[11 + 5 * j] = Math.min(r, Math.floor(x / (pj * pj)));
          bp[12 + 5 * j] = isqrtJ(xp);
        }
      }

      // large[i] for i in [1, imaxj] of block prime j, imaxj ≤ min(I0, IMAX[j]).
      // A read at i·p collides with this pass's writes only when i·p ≤ imaxj,
      // i.e. i ≤ H; ascending order is always correct, so the coordinator
      // handles the small low bands itself and the threads take over once a
      // band is big enough.  Reads above I0 are never written in this pass.
      function prefixPass(j, imaxj) {
        var p = bp[8 + 5 * j];
        if (imaxj < PAR_MIN) { prefixRange(j, 1, imaxj); return; }
        var iSw = Math.min(imaxj, Math.floor(r / p));
        var H = Math.min(iSw, Math.floor(imaxj / p));
        var seqEnd = Math.min(H, Math.max(p - 1, SEQ_MIN));
        prefixRange(j, 1, seqEnd);
        var lo = seqEnd + 1;
        while (lo <= H) {                                      // bands [lo, lo·p) ascending
          var hi = Math.min(lo * p - 1, H);
          dispatch(T_PREFIX, j, lo, hi, 0);
          lo = hi + 1;
        }
        if (H + 1 <= imaxj) {                                  // hazard-free tail
          if (imaxj - H >= DYN_MIN) dispatch(T_PREFIX, j, H + 1, imaxj, chunkFor(imaxj - H));
          else prefixRange(j, H + 1, imaxj);
        }
      }
      // small[]: descending bands, each reading only below itself; the big
      // bands (near r) go to the threads, the small tail stays sequential.
      function smallPass(p, sp1) {
        var p2 = p * p;
        if (p2 > r) return;
        var U = r;
        if (r >= PAR_MIN) {
          while (U >= p2) {
            var L = Math.max(p2, Math.floor(U / p) + 1);
            if (U - L + 1 < SEQ_MIN) break;
            dispatch(T_SMALL, L, U, p, sp1);
            U = L - 1;
          }
        }
        if (U >= p2) {
          if (kern) kern.small_range(SO, BigInt(p2), BigInt(U), BigInt(p), BigInt(sp1));
          else smallJS(p2, U, p, sp1);
        }
      }

      dispatch(T_INIT, 0, 0, 0, 0);
      // blocks start once the prefix I0 = max(r/p, x/p³) is at most r/2, i.e. p³ ≥ 2x/r
      var pb3 = 2 * Math.floor(x / r), pblock = 3;
      while (pblock * pblock * pblock < pb3) pblock++;
      var nextProg = 0, j;
      for (var p = 2; p <= r;) {
        if (p >= nextProg) { nextProg = p + 2048; post({ progress: p / r }); }
        if (small[p] === small[p - 1]) { p++; continue; }
        if (p < pblock) {
          // classical single-prime pass: a block of one with I0 = r (no sweep)
          P[0] = p;
          setBlock(1, r);
          prefixPass(0, bp[11]);
          smallPass(p, bp[10]);
          p++;
          continue;
        }
        var k = 0;
        for (var q = p; q <= r && q <= 2 * p && k < KMAX; q++) if (small[q] !== small[q - 1]) P[k++] = q;
        var p1 = P[0], I0 = Math.floor(r / p1);
        if (p1 * p1 * p1 <= x) I0 = Math.max(I0, Math.floor(x / (p1 * p1 * p1)));
        setBlock(k, I0);
        var sweepHi = bp[11];                                   // IMAX[0]
        if (sweepHi > I0) {
          if (sweepHi - I0 >= DYN_MIN) dispatch(T_SWEEP, I0 + 1, sweepHi, chunkFor(sweepHi - I0), 0);
          else if (kern) kern.sweep_range(SO, LO, BP, DO, BigInt(I0 + 1), BigInt(sweepHi));
          else sweepJS(I0 + 1, sweepHi);
        }
        if (I0 < PAR_MIN && p1 * p1 > r) {
          // every prefix is coordinator-only and there are no small-passes: one call
          if (kern) kern.block_tail(SO, LO, BP); else blockTailJS();
        } else {
          for (j = 0; j < k; j++) {
            prefixPass(j, Math.min(I0, bp[11 + 5 * j]));
            smallPass(P[j], bp[10 + 5 * j]);
          }
        }
        p = P[k - 1] + 1;
      }
      par[0] = T_EXIT;
      Atomics.add(ctrl, GEN, 1);
      Atomics.notify(ctrl, GEN);
      post({ done: true, value: large[1], kernels: kern ? (m.simd ? "wasm+simd" : "wasm") : "js" });
    });
  }

  var WORKER_SRC = "(" + workerMain.toString() + ")()";

  // the thread-safe kernels, from engine-wasm.js (inlined in the page, required in Node)
  function wasmModule() {
    try {
      var W = typeof NthPrimeWasm !== "undefined" ? NthPrimeWasm : (IS_NODE ? require("./engine-wasm.js") : null);
      return W && W.parB64 ? W : null;
    } catch (e) { return null; }
  }
  var lastKernels = "none";

  function isqrt(n) {
    var r = Math.floor(Math.sqrt(n));
    while ((r + 1) * (r + 1) <= n) r++;
    while (r * r > n) r--;
    return r;
  }

  function available() {
    if (typeof SharedArrayBuffer !== "function" || typeof Atomics !== "object") return false;
    if (IS_NODE) { try { require("worker_threads"); return true; } catch (e) { return false; } }
    return typeof Worker === "function" && typeof crossOriginIsolated !== "undefined" && crossOriginIsolated === true;
  }

  function threads() {
    var n = 1;
    if (IS_NODE) { try { n = require("os").cpus().length; } catch (e) { n = 2; } }
    else if (typeof navigator !== "undefined" && navigator.hardwareConcurrency) n = navigator.hardwareConcurrency;
    return Math.max(1, Math.min(32, n | 0));
  }

  function spawn() {
    if (IS_NODE) {
      var W = require("worker_threads").Worker;
      return { w: new W(WORKER_SRC, { eval: true }), url: null };
    }
    var url = URL.createObjectURL(new Blob([WORKER_SRC], { type: "application/javascript" }));
    return { w: new Worker(url), url: url };
  }

  // exact π(x) on K threads → Promise<number>; the promise has .cancel()
  function primeCountParallel(x, opts) {
    opts = opts || {};
    var handles = [], settled = false, rejectFn = null;
    function cleanup() {
      for (var i = 0; i < handles.length; i++) {
        try { handles[i].w.terminate(); } catch (e) { /* already gone */ }
        if (handles[i].url) { try { URL.revokeObjectURL(handles[i].url); } catch (e2) { /* ignore */ } }
      }
      handles = [];
    }
    var promise = new Promise(function (resolve, reject) {
      rejectFn = reject;
      x = Math.floor(Number(x));
      if (!(x >= 0) || x > 9e15) return reject(new RangeError("x must be an integer in [0, 9×10^15]"));
      if (x < 2) return resolve(0);
      if (!available()) return reject(new Error("multi-core engine unavailable in this environment"));
      var K = Math.max(1, Math.min(32, (opts.threads | 0) || threads()));
      if (K < 2) return reject(new Error("multi-core engine needs at least 2 threads"));
      var r = isqrt(x);
      var init = null;
      var W = opts.kernels === "js" ? null : wasmModule();
      var SEG = (W && W.SEG) || 4096;
      if (W && typeof WebAssembly === "object") {
        // one shared wasm memory holds control words, block parameters, both
        // tables and one difference array per thread; everything starts above
        // the module's 64 KiB shadow stack + globals region
        var ctrlOff = 131072, parOff = 131200, bpOff = 131328, smallOff = 133120;
        var largeOff = smallOff + 4 * (r + 2);
        largeOff += (8 - (largeOff % 8)) % 8;
        var dOff = largeOff + 8 * (r + 2);
        var total = dOff + K * 8 * (SEG + 2) + 65536;
        // exactly what we need (never grows); a 4 GiB reservation makes some
        // browsers refuse shared memory and silently fall back to JS kernels
        var pages = Math.max(17, Math.ceil(total / 65536));
        try {
          var memory = new WebAssembly.Memory({ initial: pages, maximum: pages, shared: true });
          init = { memory: memory, ctrlOff: ctrlOff, parOff: parOff, bpOff: bpOff, smallOff: smallOff, largeOff: largeOff,
                   dOff: dOff, parB64: W.parB64, simd: !!W.simd, seg: SEG };
        } catch (e) { init = null; }
      }
      if (!init) {
        try {
          init = {
            ctrl: new SharedArrayBuffer(64), par: new SharedArrayBuffer(64), bp: new SharedArrayBuffer(8 * BP_LEN),
            small: new SharedArrayBuffer(4 * (r + 2)), large: new SharedArrayBuffer(8 * (r + 2)), seg: SEG
          };
        } catch (e) {
          return reject(new Error("not enough memory for the multi-core tables (" + Math.round(12 * r / 1048576) + " MB)"));
        }
      }
      function finish(err, value) {
        if (settled) return;
        settled = true;
        cleanup();
        if (err) reject(err); else resolve(value);
      }
      for (var id = 0; id < K; id++) {
        var h;
        try { h = spawn(); } catch (e) { return finish(e); }
        handles.push(h);
        (function (w, id) {
          var onMsg = function (m) {
            if (m && m.progress !== undefined) { if (opts.onProgress && !settled) opts.onProgress(m.progress); return; }
            if (m && m.done) { lastKernels = m.kernels || "js"; finish(null, m.value); }
          };
          var onErr = function (e) {
            finish(new Error("multi-core thread " + id + " failed: " + ((e && (e.message || e.reason)) || e)));
          };
          if (IS_NODE) { w.on("message", onMsg); w.on("error", onErr); }
          else { w.onmessage = function (e) { onMsg(e.data); }; w.onerror = onErr; }
        })(h.w, id);
      }
      for (id = 0; id < K; id++) {
        var msg = { K: K, id: id, x: x, r: r };
        for (var key in init) if (Object.prototype.hasOwnProperty.call(init, key)) msg[key] = init[key];
        handles[id].w.postMessage(msg);
      }
    });
    promise.cancel = function () {
      if (settled) return;
      settled = true;
      cleanup();
      rejectFn(new Error("cancelled"));
    };
    return promise;
  }

  return {
    primeCountParallel: primeCountParallel,
    available: available,
    threads: threads,
    kernels: function () { return lastKernels; }, // "wasm+simd", "wasm" or "js" for the last completed count
    version: "2.0.0"
  };
});
