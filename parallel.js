/*
 * parallel.js — multi-core exact prime counting: the verified Lucy_Hedgehog
 * recurrence over SharedArrayBuffer, K threads, no locks in the hot loops.
 *
 * Why it is exact.  Within one prime p the recurrence writes large[i] while
 * reading large[i·p] (a HIGHER index) and writes small[v] while reading
 * small[⌊v/p⌋] (a LOWER index).  Processing the index range in geometric
 * bands — [p^k, p^{k+1}) ascending for `large`, (r/p^{k+1}, r/p^k] descending
 * for `small` — makes every band read only indices outside itself that no
 * band has written yet in this pass, so a band can be split across threads
 * freely.  One barrier per band; a few thousand barriers per query.  Primes
 * with little work (imax < 2^15) are done by the coordinator alone, exactly
 * as in the single-thread engine.  Results are cross-checked against the
 * three single-thread engines in the test suite.
 *
 * Runs in Node (worker_threads) and in browsers that are cross-origin
 * isolated (the app's service worker arranges that on GitHub Pages).  The
 * same source runs in every thread; thread 0 coordinates.
 */
(function (root, factory) {
  "use strict";
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.NthPrimeParallel = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  var IS_NODE = typeof process !== "undefined" && !!(process.versions && process.versions.node) &&
                typeof require === "function";

  // ---- thread body: self-contained, serialised with toString() ----------
  function workerMain() {
    var GEN = 0, REM = 1;
    var T_INIT = 1, T_LARGE = 2, T_SMALL = 3, T_EXIT = 9;
    var PAR_MIN = 1 << 16;   // primes with less work than this: coordinator alone
    var SEQ_MIN = 1 << 12;   // bands smaller than this: not worth waking the threads
    var isNode = typeof self === "undefined";
    var port = isNode ? require("worker_threads").parentPort : null;
    function post(m) { if (isNode) port.postMessage(m); else self.postMessage(m); }
    function onMessage(f) { if (isNode) port.on("message", f); else self.onmessage = function (e) { f(e.data); }; }

    onMessage(function (m) {
      var K = m.K, id = m.id, x = m.x, r = m.r;
      var ctrl, par, small, large, kern = null;
      if (m.memory) {
        // one shared WebAssembly memory: typed-array views for JS, offsets for the kernels
        var buf = m.memory.buffer;
        ctrl = new Int32Array(buf, m.ctrlOff, 16);
        par = new Float64Array(buf, m.parOff, 8);
        small = new Uint32Array(buf, m.smallOff, r + 2);
        large = new Float64Array(buf, m.largeOff, r + 2);
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
        small = new Uint32Array(m.small);
        large = new Float64Array(m.large);
      }
      var SO = m.smallOff || 0, LO = m.largeOff || 0, X = BigInt(x);

      // this thread's slice of the index range [lo, hi]
      function share(lo, hi) {
        var len = hi - lo + 1;
        if (len <= 0) return null;
        var per = Math.ceil(len / K);
        var a = lo + id * per;
        var b = Math.min(hi, a + per - 1);
        return a <= b ? [a, b] : null;
      }

      function task(type) {
        var c, i, v;
        if (type === T_INIT) {
          c = share(1, r);
          if (!c) return;
          if (kern) { kern.init_range(SO, LO, X, BigInt(c[0]), BigInt(c[1])); return; }
          for (v = c[0]; v <= c[1]; v++) { small[v] = v - 1; large[v] = Math.floor(x / v) - 1; }
        } else if (type === T_LARGE) {
          var p = par[3], sp1 = par[4], xp = par[5], iSw = par[6];
          c = share(par[1], par[2]);
          if (!c) return;
          if (kern) { kern.large_range(SO, LO, BigInt(c[0]), BigInt(c[1]), BigInt(p), BigInt(sp1), BigInt(xp), BigInt(iSw)); return; }
          var e1 = Math.min(c[1], iSw);
          for (i = c[0]; i <= e1; i++) large[i] -= large[i * p] - sp1;
          for (; i <= c[1]; i++) large[i] -= small[Math.floor(xp / i)] - sp1;
        } else if (type === T_SMALL) {
          var pp = par[3], s1 = par[4];
          c = share(par[1], par[2]);
          if (!c) return;
          if (kern) { kern.small_range(SO, BigInt(c[0]), BigInt(c[1]), BigInt(pp), BigInt(s1)); return; }
          var a = c[0];
          for (v = c[1]; v >= a;) {
            var q = Math.floor(v / pp);
            var w = Math.max(q * pp, a);
            var sub = small[q] - s1;
            for (; v >= w; v--) small[v] -= sub;
          }
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
      function dispatch(type, a1, a2, a3, a4, a5, a6) {
        par[0] = type; par[1] = a1; par[2] = a2; par[3] = a3; par[4] = a4; par[5] = a5; par[6] = a6;
        Atomics.store(ctrl, REM, K - 1);
        Atomics.add(ctrl, GEN, 1);
        Atomics.notify(ctrl, GEN);
        task(type);
        var rem;
        while ((rem = Atomics.load(ctrl, REM)) > 0) Atomics.wait(ctrl, REM, rem);
      }

      dispatch(T_INIT, 0, 0, 0, 0, 0, 0);
      var nextProg = 0;
      for (var p = 2; p <= r; p++) {
        if (p >= nextProg) { nextProg = p + 2048; post({ progress: p / r }); }
        if (small[p] === small[p - 1]) continue;
        var sp1 = small[p - 1], p2 = p * p, xp = Math.floor(x / p);
        var imax = Math.min(r, Math.floor(x / p2));
        var iSw = Math.min(imax, Math.floor(r / p));
        var i, v;
        if (imax < PAR_MIN) {
          // little work: identical to the single-thread engine
          if (kern) {
            kern.large_range(SO, LO, BigInt(1), BigInt(imax), BigInt(p), BigInt(sp1), BigInt(xp), BigInt(iSw));
            if (p2 <= r) kern.small_range(SO, BigInt(p2), BigInt(r), BigInt(p), BigInt(sp1));
            continue;
          }
          for (i = 1; i <= iSw; i++) large[i] -= large[i * p] - sp1;
          for (; i <= imax; i++) large[i] -= small[Math.floor(xp / i)] - sp1;
          if (p2 <= r) {
            for (v = r; v >= p2;) {
              var q = Math.floor(v / p), sub = small[q] - sp1, w = Math.max(q * p, p2);
              for (; v >= w; v--) small[v] -= sub;
            }
          }
          continue;
        }
        // large[]: a read at i·p collides with this pass's writes only when i·p ≤ imax, i.e. i ≤ H.
        // Ascending sequential order is always correct, so the coordinator handles the
        // small low bands itself and threads take over once a band is big enough.
        var H = Math.min(iSw, Math.floor(imax / p));
        var seqEnd = Math.min(H, Math.max(p - 1, SEQ_MIN));
        for (i = 1; i <= seqEnd; i++) large[i] -= large[i * p] - sp1;
        var lo = seqEnd + 1;
        while (lo <= H) {                                      // bands [lo, lo·p) ascending
          var hi = Math.min(lo * p - 1, H);
          dispatch(T_LARGE, lo, hi, p, sp1, xp, iSw);
          lo = hi + 1;
        }
        if (H + 1 <= imax) {                                   // hazard-free tail
          if (imax - H >= SEQ_MIN) dispatch(T_LARGE, H + 1, imax, p, sp1, xp, iSw);
          else {
            var t1 = Math.min(imax, iSw);
            for (i = H + 1; i <= t1; i++) large[i] -= large[i * p] - sp1;
            for (; i <= imax; i++) large[i] -= small[Math.floor(xp / i)] - sp1;
          }
        }
        // small[]: descending bands, each reading only below itself; the big
        // bands (near r) go to the threads, the small tail stays sequential.
        if (p2 <= r) {
          var U = r;
          while (U >= p2) {
            var L = Math.max(p2, Math.floor(U / p) + 1);
            if (U - L + 1 < SEQ_MIN) break;
            dispatch(T_SMALL, L, U, p, sp1, 0, 0);
            U = L - 1;
          }
          for (v = U; v >= p2;) {                              // exact single-thread order
            var q2 = Math.floor(v / p), sub2 = small[q2] - sp1, w2 = Math.max(q2 * p, p2);
            for (; v >= w2; v--) small[v] -= sub2;
          }
        }
      }
      par[0] = T_EXIT;
      Atomics.add(ctrl, GEN, 1);
      Atomics.notify(ctrl, GEN);
      post({ done: true, value: large[1], kernels: kern ? "wasm" : "js" });
    });
  }

  var WORKER_SRC = "(" + workerMain.toString() + ")()";

  // base64 of the thread-safe kernels, from engine-wasm.js (inlined in the page, required in Node)
  function parModuleB64() {
    try {
      var W = typeof NthPrimeWasm !== "undefined" ? NthPrimeWasm : (IS_NODE ? require("./engine-wasm.js") : null);
      return W && W.parB64 ? W.parB64 : "";
    } catch (e) { return ""; }
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
      var b64 = opts.kernels === "js" ? "" : parModuleB64();
      if (b64 && typeof WebAssembly === "object") {
        // one shared wasm memory holds control words + both tables
        // tables start above the module's 64 KiB shadow stack + globals region
        var ctrlOff = 131072, parOff = 131200, smallOff = 132096;
        var largeOff = smallOff + 4 * (r + 2);
        largeOff += (8 - (largeOff % 8)) % 8;
        var total = largeOff + 8 * (r + 2) + 65536;
        var pages = Math.ceil(total / 65536);
        try {
          var memory = new WebAssembly.Memory({ initial: pages, maximum: 65535, shared: true });
          init = { memory: memory, ctrlOff: ctrlOff, parOff: parOff, smallOff: smallOff, largeOff: largeOff, parB64: b64 };
        } catch (e) { init = null; }
      }
      if (!init) {
        try {
          init = {
            ctrl: new SharedArrayBuffer(64), par: new SharedArrayBuffer(64),
            small: new SharedArrayBuffer(4 * (r + 2)), large: new SharedArrayBuffer(8 * (r + 2))
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
    kernels: function () { return lastKernels; }, // "wasm" or "js" for the last completed count
    version: "1.1.0"
  };
});
