/*
 * nthprime.js — compute the n-th prime number, fast and exactly,
 * for every 1 ≤ n ≤ 2×10^14 (answers approach the 2^53 float64 frontier).
 *
 * TWO independent exact prime-counting engines live in this file:
 *   · Lucy_Hedgehog tables (O(x^{3/4}), the speed engine — primeCount)
 *   · Lagarias–Miller–Odlyzko (O(x^{2/3})-class, the verification engine —
 *     primeCountLMO: ordinary/special Möbius leaves + Fenwick segmented
 *     sieve + P₂ sweep).  Two unrelated algorithms, two unrelated bug
 *     surfaces; the test suite requires digit-for-digit agreement.
 *
 * Strategy (the same one used by state-of-the-art tools such as primecount):
 *
 *   1. ESTIMATE  x0 ≈ p(n) by inverting the Riemann prime-counting
 *      approximation R(x) with Newton's method.  R(x) is evaluated with the
 *      Gram series  R(x) = 1 + Σ_{k≥1} (ln x)^k / (k·k!·ζ(k+1)),  and the
 *      Newton iteration is seeded with Cipolla's 1902 asymptotic expansion
 *      of p(n).  Empirically |π(x0) − n| is tiny (≈ a few thousand even at
 *      n = 10^12), and the estimate is clamped into the rigorous
 *      Rosser/Dusart bracket  n(ln n + ln ln n − 1) ≤ p(n) ≤ n(ln n + ln ln n)
 *      so correctness never depends on the quality of the guess.
 *
 *   2. COUNT  π(x0) *exactly* with Lucy_Hedgehog's algorithm — a
 *      Meissel-style dynamic programme over the O(√x) distinct values of
 *      ⌊x/k⌋ that runs in O(x^{3/4}) time and O(√x) space.  No factoring,
 *      no primality tests, just arithmetic.
 *
 *   3. WALK  the remaining distance (n − π(x0) primes) with an odds-only,
 *      segmented sieve of Eratosthenes.  Because step 1 is so accurate this
 *      walk covers only a few hundred thousand integers at worst.
 *
 * Everything stays inside IEEE-754 double exact-integer range (< 2^53).
 * Supported domain: 1 ≤ n ≤ 10^12   (p(10^12) = 29 996 224 275 833 ≈ 3·10^13).
 *
 * Exactness of Math.floor(a / b) for integers: IEEE-754 division is
 * correctly rounded, so the computed quotient differs from the true rational
 * a/b by at most ½ulp ≤ (a/b)·2^-53.  A non-integer rational a/b is at least
 * 1/b away from the nearest integer, and (a/b)·2^-53 < 1/b  ⇔  a < 2^53.
 * Hence floor/ceil of a quotient of exact integer doubles below 2^53 is
 * always exact.  idiv() below still guards the outer-loop divisions out of
 * an abundance of caution (the guard is a never-taken branch in practice).
 *
 * Works in browsers (window.NthPrime), Web Workers (self.NthPrime) and
 * Node.js (module.exports).  Zero dependencies.
 */
(function (root, factory) {
  "use strict";
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.NthPrime = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  var MAX_N = 2e14; // p(MAX_N) ≈ 7.34e15, still below the 2^53 exact-integer bound
  var SIEVE_PATH_MAX = 500000; // below this a direct sieve is instant
  var FIRST_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

  // ------------------------------------------------------------------
  // exact integer helpers
  // ------------------------------------------------------------------

  function isqrt(n) {
    var r = Math.floor(Math.sqrt(n));
    while ((r + 1) * (r + 1) <= n) r++;
    while (r * r > n) r--;
    return r;
  }

  // floor(a/b) with a belt-and-suspenders correction step.
  function idiv(a, b) {
    var q = Math.floor(a / b);
    if (q * b > a) q--;
    else if ((q + 1) * b <= a) q++;
    return q;
  }

  // ------------------------------------------------------------------
  // simple sieve: all primes ≤ limit (odds-only Eratosthenes)
  // ------------------------------------------------------------------

  function primesUpTo(limit) {
    limit = Math.floor(limit);
    if (limit < 2) return new Int32Array(0);
    if (limit === 2) return Int32Array.from([2]);
    var m = (limit - 1) >>> 1; // odd number 2i+1 ↔ index i, 1 ≤ i ≤ m
    var comp = new Uint8Array(m + 1);
    var sq = isqrt(limit);
    for (var i = 1; 2 * i + 1 <= sq; i++) {
      if (!comp[i]) {
        var p = 2 * i + 1;
        for (var j = (p * p - 1) >>> 1; j <= m; j += p) comp[j] = 1;
      }
    }
    var cnt = 1;
    for (i = 1; i <= m; i++) if (!comp[i]) cnt++;
    var out = new Int32Array(cnt);
    out[0] = 2;
    var k = 1;
    for (i = 1; i <= m; i++) if (!comp[i]) out[k++] = 2 * i + 1;
    return out;
  }

  // Rosser's theorem: p(n) < n(ln n + ln ln n) for n ≥ 6.
  function upperBoundForNthPrime(n) {
    if (n < 6) return 13;
    var L = Math.log(n);
    return Math.ceil(n * (L + Math.log(L))) + 10;
  }

  // Dusart (1999): p(n) ≥ n(ln n + ln ln n − 1) for n ≥ 2.
  function lowerBoundForNthPrime(n) {
    if (n < 2) return 2;
    var L = Math.log(n);
    return Math.max(2, Math.floor(n * (L + Math.log(L) - 1)));
  }

  function nthPrimeBySieve(n) {
    return primesUpTo(upperBoundForNthPrime(n))[n - 1];
  }

  // ------------------------------------------------------------------
  // Riemann R function via the Gram series, and its inverse
  // ------------------------------------------------------------------

  // ζ(s) for real s ≥ 2 by Euler–Maclaurin; absolute error ≲ 1e-15.
  function zeta(s) {
    if (s > 64) return 1 + Math.pow(2, -s) + Math.pow(3, -s);
    var J = 24;
    var sum = 0;
    for (var j = 1; j <= J; j++) sum += Math.pow(j, -s);
    sum += Math.pow(J, 1 - s) / (s - 1) - 0.5 * Math.pow(J, -s);
    sum += (s / 12) * Math.pow(J, -s - 1);
    sum -= ((s * (s + 1) * (s + 2)) / 720) * Math.pow(J, -s - 3);
    sum += ((s * (s + 1) * (s + 2) * (s + 3) * (s + 4)) / 30240) * Math.pow(J, -s - 5);
    return sum;
  }

  var zetaCache = [];
  function zetaInt(m) {
    var z = zetaCache[m];
    if (z === undefined) {
      z = zeta(m);
      zetaCache[m] = z;
    }
    return z;
  }

  // R(x) = 1 + Σ_{k≥1} (ln x)^k / (k · k! · ζ(k+1))   (Gram series)
  function riemannR(x) {
    if (!(x > 1)) return 0;
    var L = Math.log(x);
    var sum = 1;
    var u = 1; // u = L^k / k!
    for (var k = 1; k < 500; k++) {
      u *= L / k;
      var term = u / (k * zetaInt(k + 1));
      sum += term;
      if (k > L && term < sum * 1e-16) break;
    }
    return sum;
  }

  // Cipolla (1902): p(n) ≈ n·(L + LL − 1 + (LL−2)/L − (LL² − 6LL + 11)/(2L²))
  // with L = ln n, LL = ln ln n.
  function nthPrimeEstimate(n) {
    if (n < 6) return [2, 2, 3, 5, 7, 11][n] || 2;
    var L = Math.log(n);
    var LL = Math.log(L);
    return (
      n * (L + LL - 1 + (LL - 2) / L - (LL * LL - 6 * LL + 11) / (2 * L * L))
    );
  }

  // Solve R(x) = n for x with Newton's method (R'(x) ≈ 1/ln x).
  function inverseRiemannR(n) {
    var x = nthPrimeEstimate(n);
    if (n < 6) return x;
    for (var i = 0; i < 60; i++) {
      var dx = (riemannR(x) - n) * Math.log(x);
      x -= dx;
      if (Math.abs(dx) < 0.5) break;
    }
    return x;
  }

  // Round the R-inverse guess and clamp it into the rigorous bracket.
  function guessNthPrime(n) {
    var x = Math.round(inverseRiemannR(n));
    var lo = lowerBoundForNthPrime(n);
    var hi = upperBoundForNthPrime(n);
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return x;
  }

  // ------------------------------------------------------------------
  // exact prime counting — Lucy_Hedgehog, O(x^{3/4}) time, O(√x) space
  // ------------------------------------------------------------------
  //
  // S(v) starts as |[2..v]| = v−1 and after processing prime p holds the
  // count of integers in [2..v] that are prime OR have all prime factors
  // greater than p.  Sieving by p removes exactly  S(v/p) − S(p−1)  values
  // (composites whose least prime factor is p), giving the recurrence
  //     S(v) ← S(v) − ( S(⌊v/p⌋) − π(p−1) )        for v ≥ p².
  // Only the O(√x) distinct values v = ⌊x/k⌋ ever matter:  small[v] stores
  // S(v) for v ≤ √x and large[i] stores S(⌊x/i⌋) for i ≤ √x.  After all
  // p ≤ √x are processed, large[1] = S(x) = π(x).

  // ------------------------------------------------------------------
  // exact prime counting — Lagarias–Miller–Odlyzko, ~O(x^{2/3}) time
  // ------------------------------------------------------------------
  //
  // The Meissel identity with a = π(y), x^{1/3} ≤ y < √x:
  //
  //     π(x) = φ(x, a) + a − 1 − P₂(x, a)
  //
  //   · φ(x, a)  = #{m ≤ x : no prime factor ≤ p_a}   (counts 1)
  //   · P₂(x, a) = #{m ≤ x : m = p·q, primes y < p ≤ q}
  //              = Σ_{y < p ≤ √x} ( π(x/p) − π(p) + 1 )
  //
  // P₂ needs π(v) only at v = ⌊x/p⌋ ≤ x/y: one segmented-sieve sweep up
  // to x/y answers every query in ascending order.
  //
  // φ(x, a) follows LMO: expand φ(u, b) = φ(u, b−1) − φ(u/p_b, b−1)
  // into a tree.  A node is (n, b) with term μ(n)·φ(x/n, b), n squarefree,
  // built from distinct primes of index > b.  Children of (n, b) are
  // (n, b−1) and (n·p_b, b−1).  Stop at:
  //   · ordinary leaves  (n, c):     μ(n)·φ(x/n, c) with n ≤ y — evaluated
  //     with a wheel mod W = p_1…p_c (we use c = 7, W = 510510);
  //   · special leaves   n·p_b > y:  −μ(n)·φ(x/(n·p_b), b−1).
  // Every special-leaf argument v = ⌊x/(n·p_b)⌋ is < x/y, so a second
  // segmented sieve over [1, x/y] — processing primes p_{c+1}…p_a in order
  // and answering "count of unmarked ≤ v" via a Fenwick tree right before
  // sieving p_b — evaluates all of them.
  //
  // y = α·x^{1/3} trades leaf count against sieve length x/y (primecount's
  // α tuning); correctness holds for any x^{1/3} ≤ y < √x.
  //
  // Everything is verified at runtime against the independent O(x^{3/4})
  // engine above on overlapping ranges (see test/test.js).

  function icbrt(n) {
    var r = Math.floor(Math.cbrt(n));
    while ((r + 1) * (r + 1) * (r + 1) <= n) r++;
    while (r * r * r > n) r--;
    return r;
  }

  // wheel for φ(u, 7): W = 2·3·5·7·11·13·17, tab[i] = #{1 ≤ j ≤ i, gcd(j,W)=1}
  var WHEEL_C = 7;
  var WHEEL_PRIMES = [2, 3, 5, 7, 11, 13, 17];
  var wheelCache = null;
  function getWheel() {
    if (wheelCache) return wheelCache;
    var W = 510510;
    var coprime = new Uint8Array(W + 1);
    coprime.fill(1);
    coprime[0] = 0;
    for (var k = 0; k < WHEEL_PRIMES.length; k++) {
      for (var m = WHEEL_PRIMES[k]; m <= W; m += WHEEL_PRIMES[k]) coprime[m] = 0;
    }
    var tab = new Int32Array(W + 1);
    var c = 0;
    for (var i = 1; i <= W; i++) {
      c += coprime[i];
      tab[i] = c;
    }
    wheelCache = { W: W, phiW: c, tab: tab };
    return wheelCache;
  }

  function phi7(u) {
    // φ(u, 7) — exact for 0 ≤ u < 2^53 (q·phiW ≤ u stays exact)
    if (u <= 0) return 0;
    var w = getWheel();
    var q = idiv(u, w.W);
    return q * w.phiW + w.tab[u - q * w.W];
  }

  // μ(n) and least-prime-factor for n ≤ y by sieve
  function muLpfTables(y) {
    var lpf = new Int32Array(y + 1);
    for (var i = 2; i <= y; i++) {
      if (lpf[i] === 0) {
        for (var j = i; j <= y; j += i) if (lpf[j] === 0) lpf[j] = i;
      }
    }
    var mu = new Int8Array(y + 1);
    mu.fill(1);
    mu[0] = 0;
    for (var p = 2; p <= y; p++) {
      if (lpf[p] === p) {
        for (var m = p; m <= y; m += p) mu[m] = -mu[m];
        var pp = p * p;
        if (pp <= y) for (var m2 = pp; m2 <= y; m2 += pp) mu[m2] = 0;
      }
    }
    return { mu: mu, lpf: lpf };
  }

  // Fenwick tree over [1..size] supporting point add / prefix sum
  function Fenwick(size) {
    this.n = size;
    this.t = new Int32Array(size + 1);
  }
  Fenwick.prototype.buildFromUnit = function (unit, len) {
    // unit: Uint8Array (1 = present), positions 1..len
    var t = this.t;
    t.fill(0);
    var n = this.n;
    for (var i = 1; i <= len; i++) {
      t[i] += unit[i];
      var j = i + (i & -i);
      if (j <= n) t[j] += t[i];
    }
  };
  Fenwick.prototype.remove = function (i) {
    for (; i <= this.n; i += i & -i) this.t[i]--;
  };
  Fenwick.prototype.prefix = function (i) {
    var s = 0;
    for (; i > 0; i -= i & -i) s += this.t[i];
    return s;
  };

  // radix sort of Uint32 keys (returns order as Uint32Array of indices)
  function sortIndicesByKey(keys) {
    var n = keys.length;
    var order = new Uint32Array(n);
    var tmp = new Uint32Array(n);
    for (var i = 0; i < n; i++) order[i] = i;
    var count = new Uint32Array(65536);
    for (var shift = 0; shift <= 16; shift += 16) {
      count.fill(0);
      for (i = 0; i < n; i++) count[(keys[order[i]] >>> shift) & 65535]++;
      var pos = 0;
      for (i = 0; i < 65536; i++) {
        var c = count[i];
        count[i] = pos;
        pos += c;
      }
      for (i = 0; i < n; i++) tmp[count[(keys[order[i]] >>> shift) & 65535]++] = order[i];
      var sw = order;
      order = tmp;
      tmp = sw;
    }
    return order;
  }

  // P₂(x, a) = Σ_{y < p ≤ √x} (π(x/p) − π(p) + 1), one ascending sweep.
  // basePrimes must contain every prime ≤ √x.
  function p2Count(x, y, basePrimes, onProgress) {
    var sqrtX = isqrt(x);
    // prime indices: π(basePrimes[k]) = k + 1
    var k1 = basePrimes.length - 1;
    while (k1 >= 0 && basePrimes[k1] > sqrtX) k1--;
    var k0 = k1;
    while (k0 >= 0 && basePrimes[k0] > y) k0--;
    k0++; // first index with p > y
    if (k0 > k1) return 0;

    var total = 0;
    var cum = 1; // π so far in the sweep (prime 2 precounted)
    var k = k1;  // descending p ⇒ ascending v = ⌊x/p⌋
    var vmax = idiv(x, basePrimes[k0]);
    var SEG = 1 << 21;
    var lo = 3;
    var done = 0;
    while (lo <= vmax && k >= k0) {
      var hi = Math.min(lo + SEG - 1, vmax);
      var seg = sieveWindow(lo, hi, basePrimes);
      var comp = seg.comp, base = seg.base, len = seg.len;
      var idx = 0;
      while (k >= k0) {
        var v = idiv(x, basePrimes[k]);
        if (v > hi) break;
        var target = v < base ? -1 : (v - base) >> 1; // floor → last odd index ≤ v
        while (idx <= target) {
          if (comp[idx] === 0) cum++;
          idx++;
        }
        // π(v) = cum; π(p) = k + 1; term = π(v) − π(p) + 1
        total += cum - (k + 1) + 1;
        k--;
      }
      while (idx < len) {
        if (comp[idx] === 0) cum++;
        idx++;
      }
      done += hi - lo + 1;
      if (onProgress) onProgress(done / vmax);
      lo = hi + 1;
    }
    if (k >= k0) throw new Error("internal error: P2 sweep ended early");
    return total;
  }

  // φ(x, a) via LMO ordinary + special leaves.  y defines a = π(y) = aCount.
  // basePrimes must cover every prime ≤ y.  1-based prime indexing:
  // π(basePrimes[k]) = k + 1, so p_b = basePrimes[b − 1].
  function phiLMO(x, y, basePrimes, aCount, onProgress) {
    var tables = muLpfTables(y);
    var mu = tables.mu, lpf = tables.lpf;
    var P_C = 17; // p_c for c = WHEEL_C = 7
    var a = aCount;

    // ---- ordinary leaves: Σ_{n ≤ y, μ(n)≠0, lpf(n) > 17 (or n=1)} μ(n)·φ(x/n, c)
    var ordinary = 0;
    for (var n = 1; n <= y; n++) {
      if (n === 1 || (mu[n] !== 0 && lpf[n] > P_C)) {
        ordinary += mu[n] * phi7(idiv(x, n));
      }
    }

    // ---- special leaves: pairs (n, p_b), μ(n)≠0, 17 < p_b < lpf(n),
    //      p_b ≤ p_a, n ≤ y < n·p_b.  Term: −μ(n)·φ(⌊x/(n·p_b)⌋, b−1).
    // n = 1 never qualifies (1·p_b ≤ p_a ≤ y).
    // Two passes: count, then fill exact-sized typed arrays.
    var loBound = new Int32Array(y + 1);
    var hiBound = new Int32Array(y + 1);
    var Q = 0;
    for (n = 2; n <= y; n++) {
      loBound[n] = 1;
      hiBound[n] = 0;
      if (mu[n] === 0) continue;
      var lp = lpf[n];
      if (lp <= P_C) continue;
      // p_b > y/n  ⟺  n·p_b > y  (integers: p_b ≥ ⌊y/n⌋+1 ⇒ n·p_b > y)
      var loK = lowerBoundPrime(basePrimes, Math.max(P_C, idiv(y, n)) + 0.5);
      var hiK = lowerBoundPrime(basePrimes, lp - 0.5) - 1; // last prime < lpf(n)
      if (hiK > a - 1) hiK = a - 1;                        // p_b ≤ p_a
      if (hiK >= loK) {
        loBound[n] = loK;
        hiBound[n] = hiK;
        Q += hiK - loK + 1;
      }
    }
    var qV = new Float64Array(Q);   // v ≤ x/y can exceed 2^32
    var qB = new Int32Array(Q);
    var qSign = new Int8Array(Q);
    var qi = 0;
    for (n = 2; n <= y; n++) {
      for (var kk = loBound[n]; kk <= hiBound[n]; kk++) {
        qV[qi] = idiv(x, n * basePrimes[kk]);
        qB[qi] = kk + 1; // 1-based index b
        qSign[qi] = -mu[n];
        qi++;
      }
    }
    var vmax = 1;
    for (var qi = 0; qi < Q; qi++) if (qV[qi] > vmax) vmax = qV[qi];

    // ---- segmented Fenwick sieve over [1, vmax] answers all queries.
    // Queries grouped by (segment, b) ascending via one radix sort:
    //   key = segId·(a+2) + b   (fits 32 bits for every supported x)
    var SEG = 1 << 20;
    var nb = a + 2;
    var segCount = Math.floor((vmax - 1) / SEG) + 1;
    if (segCount * nb > 4294967295) throw new RangeError("phiLMO: key overflow");
    var keys = new Uint32Array(Q);
    for (qi = 0; qi < Q; qi++) {
      keys[qi] = Math.floor((qV[qi] - 1) / SEG) * nb + qB[qi];
    }
    var order = Q > 0 ? sortIndicesByKey(keys) : new Uint32Array(0);

    // cnt[b] = #unmarked in [1, current segment start) at state b   (b ∈ [c, a])
    var cnt = new Float64Array(a + 1);
    // nxt[k] = next unsieved multiple of basePrimes[k] (cursor across segments)
    var firstB = WHEEL_C; // 0-based index of p_{c+1} = basePrimes[7] = 19
    var nxt = new Float64Array(a);
    for (var k2 = firstB; k2 < a; k2++) nxt[k2] = basePrimes[k2];

    var special = 0;
    var mark = new Uint8Array(SEG + 1);
    var unit = new Uint8Array(SEG + 1);
    var bit = new Fenwick(SEG);
    var qPtr = 0;

    for (var segId = 0; segId < segCount; segId++) {
      var lo = segId * SEG + 1;
      var hi = Math.min(lo + SEG - 1, vmax);
      var len = hi - lo + 1;

      // state c: mark multiples of the wheel primes (number 1 stays unmarked)
      mark.fill(0);
      for (var wk = 0; wk < WHEEL_PRIMES.length; wk++) {
        var wp = WHEEL_PRIMES[wk];
        for (var m = Math.ceil(lo / wp) * wp; m <= hi; m += wp) mark[m - lo + 1] = 1;
      }
      for (var ii = 1; ii <= len; ii++) unit[ii] = 1 - mark[ii];
      if (len < SEG) unit.fill(0, len + 1);
      bit.buildFromUnit(unit, SEG);
      var curUnm = bit.prefix(len); // unmarked in this segment at current state

      // walk states c → a: answer b-queries at state b−1, then sieve p_b
      for (var bIdx = firstB; bIdx < a; bIdx++) {
        var b = bIdx + 1;          // we are at state b−1; about to sieve p_b
        while (qPtr < Q) {
          var qq = order[qPtr];
          if (qB[qq] !== b || Math.floor((qV[qq] - 1) / SEG) !== segId) break;
          special += qSign[qq] * (cnt[b - 1] + bit.prefix(qV[qq] - lo + 1));
          qPtr++;
        }
        cnt[b - 1] += curUnm;      // state b−1 finalized for this segment
        var pb = basePrimes[bIdx];
        var mm = nxt[bIdx];
        if (mm <= hi) {
          for (; mm <= hi; mm += pb) {
            var pos = mm - lo + 1;
            if (!mark[pos]) {
              mark[pos] = 1;
              bit.remove(pos);
              curUnm--;
            }
          }
          nxt[bIdx] = mm;
        }
      }
      cnt[a] += curUnm;
      if (onProgress && (segId & 15) === 0) onProgress((segId + 1) / segCount);
    }
    if (qPtr !== Q) throw new Error("internal error: unanswered special leaves");

    return ordinary + special;
  }

  // first index k with basePrimes[k] > bound (binary search)
  function lowerBoundPrime(primes, bound) {
    var lo = 0, hi = primes.length;
    while (lo < hi) {
      var mid = (lo + hi) >>> 1;
      if (primes[mid] > bound) hi = mid;
      else lo = mid + 1;
    }
    return lo;
  }

  // y = α·x^{1/3}: bigger α shortens both sieves (range x/y) but multiplies
  // the special-leaf count (≈ π(y)²/2 pairs).  Without the Deléglise–Rivat
  // easy/trivial-leaf split, α = 1 is the sweet spot; correctness holds for
  // any α (the test suite cross-checks several).
  function defaultAlpha() {
    return 1;
  }

  // exact π(x) via the LMO machinery above
  function primeCountLMO(x, onProgress, opts) {
    x = Math.floor(x);
    if (x < 2) return 0;
    if (x > 9e15) throw new RangeError("primeCountLMO: x exceeds 2^53 safety bound");
    var sqrtX = isqrt(x);
    var basePrimes = (opts && opts.basePrimes) || primesUpTo(sqrtX);
    if (basePrimes.length < 9) return primeCount(x); // tiny x: wheel assumptions fail
    var y0 = icbrt(x);
    var alpha = (opts && opts.alpha) || defaultAlpha();
    var y = Math.min(y0 * alpha, Math.max(y0, sqrtX - 1));
    var a = lowerBoundPrime(basePrimes, y + 0.5); // π(y)
    if (a <= WHEEL_C + 1) return primeCount(x);
    var phi = phiLMO(x, y, basePrimes, a,
      onProgress ? function (f) { onProgress(f * 0.6); } : null);
    var p2 = p2Count(x, y, basePrimes,
      onProgress ? function (f) { onProgress(0.6 + f * 0.4); } : null);
    return phi + a - 1 - p2;
  }

  // ------------------------------------------------------------------
  // optional compiled engine (engine-wasm.js, generated from engine.c):
  // the same verified Lucy recurrence, ~1.7× faster, near-native speed
  // ------------------------------------------------------------------

  var wasmMod = null, wasmTried = false;
  function getWasmModule() {
    if (wasmTried) return wasmMod;
    wasmTried = true;
    try {
      var W = null;
      if (typeof NthPrimeWasm !== "undefined") W = NthPrimeWasm; // inlined page / worker
      else if (typeof require === "function") W = require("./engine-wasm.js"); // Node
      if (W && W.init && W.init()) wasmMod = W;
    } catch (e) {
      wasmMod = null;
    }
    return wasmMod;
  }

  // exact π(x) on the compiled engine; null ⇒ unavailable (caller falls back)
  function primeCountWasm(x, onProgress) {
    var W = getWasmModule();
    if (!W) return null;
    var w = W.init();
    var r = isqrt(x);
    var smallOff = 131072; // above the module's shadow stack + globals region
    var largeOff = smallOff + 4 * (r + 2);
    largeOff += (8 - (largeOff % 8)) % 8; // u64 table must be 8-aligned
    var scratchOff = largeOff + 8 * (r + 2); // the sweep's difference array (SEG + 2 u64)
    var need = scratchOff + 8 * ((W.SEG || 4096) + 2) + 65536;
    var have = w.memory.buffer.byteLength;
    if (need > have) {
      try { w.memory.grow(Math.ceil((need - have) / 65536)); }
      catch (e) { return null; } // out of memory — JS engine takes over
    }
    W.setProgress(onProgress || null);
    var v = Number(w.exports.pi_lucy(BigInt(x), smallOff, largeOff, scratchOff));
    W.setProgress(null);
    return v;
  }

  // Engine selector.  Measured head-to-head in JS, the Lucy_Hedgehog tables
  // beat LMO(α=1) by ~2× at every size we support (LMO's asymptotic edge of
  // x^{1/12} cannot overcome its Fenwick/segment constants below ~10^17),
  // so Lucy computes and LMO serves as the independent cross-checking
  // engine (see test/test.js and the `engine` option of the CLI).
  function primeCountAuto(x, onProgress, opts) {
    var eng = opts && opts.engine;
    if (eng === "lmo") return primeCountLMO(x, onProgress, opts);
    if (eng === "lucy" || eng === "js") return primeCount(x, onProgress);
    if (eng === "wasm" || x >= 1e7) {
      var v = primeCountWasm(x, onProgress);
      if (v !== null) return v;
      if (eng === "wasm") throw new Error("WebAssembly engine unavailable in this environment");
    }
    return primeCount(x, onProgress);
  }

  function primeCount(N, onProgress) {
    N = Math.floor(N);
    if (N < 2) return 0;
    if (N > 9e15) throw new RangeError("primeCount: x exceeds 2^53 safety bound");
    var r = isqrt(N);
    var small = new Uint32Array(r + 1); // S(v) ≤ v ≤ r < 2^32 — half the bytes of the hottest array
    var large = new Float64Array(r + 1); // large[i] = S(⌊N/i⌋), values up to N
    var v, i;
    for (v = 1; v <= r; v++) small[v] = v - 1;
    for (i = 1; i <= r; i++) large[i] = idiv(N, i) - 1;

    for (var p = 2; p <= r; p++) {
      if (small[p] === small[p - 1]) continue; // p is composite
      var sp1 = small[p - 1]; // π(p−1)
      var p2 = p * p;
      var Np = idiv(N, p);
      var imax = Math.min(r, idiv(N, p2));
      // ⌊N/(i·p)⌋ is large[i·p] while i·p ≤ r, else small[⌊Np/i⌋]
      // (⌊⌊N/p⌋/i⌋ = ⌊N/(p·i)⌋).
      var iSwitch = Math.min(imax, Math.floor(r / p));
      for (i = 1; i <= iSwitch; i++) large[i] -= large[i * p] - sp1;
      // beyond √Np consecutive i share one quotient q; the run of q is
      // (⌊Np/(q+1)⌋, ⌊Np/q⌋] (i ≤ ⌊Np/q⌋ ⇒ Np/i ≥ q and i > ⌊Np/(q+1)⌋ ⇒
      // Np/i < q+1), so walking q downward costs one independent division
      // per run instead of a latency chain of two.
      var iGroup = Math.min(imax + 1, Math.max(i, isqrt(Np)));
      for (; i < iGroup; i++) large[i] -= small[Math.floor(Np / i)] - sp1;
      if (i <= imax) {
        var q0 = Math.floor(Np / i), qMin = Math.max(1, Math.floor(Np / imax));
        var ePrev = i - 1;
        for (; q0 >= qMin; q0--) {
          var iEnd = Math.min(imax, Math.floor(Np / q0));
          var c0 = small[q0] - sp1;
          for (var j = ePrev + 1; j <= iEnd; j++) large[j] -= c0;
          ePrev = iEnd;
        }
      }
      // Descending update keeps small[q] (q = ⌊v/p⌋ < v) at its
      // previous-pass value until v itself reaches q.  Values of v sharing
      // the same quotient q are processed as one block.
      for (v = r; v >= p2; ) {
        var q = (v / p) | 0; // v ≤ r ≤ 5.6e6 → 32-bit safe
        var sub = small[q] - sp1;
        var w = q * p;
        if (w < p2) w = p2;
        for (; v >= w; v--) small[v] -= sub;
      }
      if (onProgress && (p & 2047) === 0) onProgress(p / r);
    }
    return large[1];
  }

  // ------------------------------------------------------------------
  // segmented sieve of a window [lo, hi] (odd numbers only)
  // ------------------------------------------------------------------
  //
  // Returns { base, len, comp } where index i represents the odd number
  // base + 2i and comp[i] === 0 means that number is prime.
  // Requires basePrimes to cover every prime ≤ √hi, and lo ≥ 3.
  // NOTE: lo/hi may exceed 2^31, so no bitwise ops on them — arithmetic only.

  function sieveWindow(lo, hi, basePrimes) {
    if (lo % 2 === 0) lo += 1;
    if (hi % 2 === 0) hi -= 1;
    if (hi < lo) return { base: lo, len: 0, comp: new Uint8Array(0) };
    var len = (hi - lo) / 2 + 1;
    var comp = new Uint8Array(len);
    var sq = isqrt(hi);
    for (var k = 1; k < basePrimes.length; k++) { // skip basePrimes[0] = 2
      var p = basePrimes[k];
      if (p > sq) break;
      var start = p * p;
      if (start < lo) {
        start = Math.ceil(lo / p) * p; // exact: lo < 2^53
        if (start % 2 === 0) start += p;
      }
      for (var j = (start - lo) / 2; j < len; j += p) comp[j] = 1;
    }
    return { base: lo, len: len, comp: comp };
  }

  function countZeros(comp, len) {
    var c = 0;
    for (var i = 0; i < len; i++) if (comp[i] === 0) c++;
    return c;
  }

  // ------------------------------------------------------------------
  // deterministic POLYNOMIAL-TIME primality test (Miller–Rabin, BigInt)
  // ------------------------------------------------------------------
  //
  // Testing whether m is prime runs in O(d³) bit operations where d is the
  // number of DIGITS of m — genuinely polynomial time in the input size
  // (the same complexity class as AKS, with far better constants).  It is
  // deterministic — a proof, not a probability — for every
  // m < 3,317,044,064,679,887,385,961,981 ≈ 3.3×10²⁴ using the first 13
  // prime witnesses (Sorenson & Webster, Math. Comp. 85, 2016).
  //
  // Note the contrast documented in the README: FINDING the n-th prime in
  // time polynomial in the digits of n is an open problem in mathematics;
  // TESTING a given number is polynomial and implemented right here.

  var BI0, BI1, BI2, MR_LIMIT, MR_WITNESSES;
  function initBigints() {
    if (BI0 !== undefined) return true;
    if (typeof BigInt === "undefined") return false;
    BI0 = BigInt(0);
    BI1 = BigInt(1);
    BI2 = BigInt(2);
    MR_LIMIT = BigInt("3317044064679887385961981");
    MR_WITNESSES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41].map(BigInt);
    return true;
  }

  function modpow(b, e, m) {
    var r = BI1;
    b %= m;
    while (e > BI0) {
      if (e % BI2 === BI1) r = (r * b) % m;
      b = (b * b) % m;
      e /= BI2;
    }
    return r;
  }

  var TRIAL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];

  // isPrime(m) → { prime: bool, factor?: string }  for 0 ≤ m < 3.3×10²⁴.
  // Accepts Number, BigInt, or a numeric string of any length in range.
  function isPrime(input) {
    if (!initBigints()) throw new Error("BigInt unavailable in this environment");
    var m = typeof input === "bigint" ? input : BigInt(String(input).trim());
    if (m < BI2) return { prime: false, reason: "primes start at 2" };
    if (m >= MR_LIMIT) throw new RangeError("primality test is deterministic only below 3.3×10^24");
    var i, p;
    for (i = 0; i < TRIAL_PRIMES.length; i++) {
      p = BigInt(TRIAL_PRIMES[i]);
      if (m === p) return { prime: true };
      if (m % p === BI0) return { prime: false, factor: String(TRIAL_PRIMES[i]) };
    }
    // m odd, no factor ≤ 97 — write m−1 = d·2^s and run the witness set
    var d = m - BI1;
    var s = 0;
    while (d % BI2 === BI0) { d /= BI2; s++; }
    for (i = 0; i < MR_WITNESSES.length; i++) {
      var a = MR_WITNESSES[i] % m;
      if (a === BI0) continue;
      var x = modpow(a, d, m);
      if (x === BI1 || x === m - BI1) continue;
      var composite = true;
      for (var r = 1; r < s; r++) {
        x = (x * x) % m;
        if (x === m - BI1) { composite = false; break; }
      }
      if (composite) return { prime: false, witness: String(MR_WITNESSES[i]) };
    }
    return { prime: true };
  }

  // ------------------------------------------------------------------
  // neighbouring primes — lets every answer carry an on-the-spot local
  // verification (p is prime; prev/next primes and gaps), a few dozen
  // deterministic Miller–Rabin tests: microseconds.
  // ------------------------------------------------------------------

  function toNumberIfSafe(b) {
    return b <= BigInt("9007199254740991") ? Number(b) : b.toString();
  }

  function primeNeighbors(p) {
    if (!initBigints()) throw new Error("BigInt unavailable in this environment");
    var P = typeof p === "bigint" ? p : BigInt(String(p));
    if (P < BI2) throw new RangeError("primeNeighbors: p must be ≥ 2");
    var q = P - BI1, prev = null;
    while (q >= BI2) {
      if (isPrime(q).prime) { prev = q; break; }
      q -= BI1;
    }
    var r = P + BI1;
    while (!isPrime(r).prime) r += BI1;
    return { prev: prev === null ? null : toNumberIfSafe(prev), next: toNumberIfSafe(r) };
  }

  function withNeighbors(res) {
    try {
      var nb = primeNeighbors(res.value);
      res.prev = nb.prev;
      res.next = nb.next;
    } catch (e) { /* neighbours are a courtesy; never fail the answer over them */ }
    return res;
  }

  // ------------------------------------------------------------------
  // π(x) as a first-class query (the page's pi(x) mode and the CLI)
  // ------------------------------------------------------------------

  var MAX_PI_X = 9e15; // exact-integer safety: below 2^53

  function countPrimes(x, opts) {
    x = typeof x === "string" ? Number(x) : x;
    if (typeof x !== "number" || !Number.isFinite(x) || !Number.isInteger(x)) {
      throw new RangeError("x must be an integer");
    }
    if (x < 0) throw new RangeError("x must be ≥ 0");
    if (x > MAX_PI_X) throw new RangeError("x must be ≤ 9×10^15 (answers stay below 2^53)");
    var onProgress = opts && opts.onProgress;
    var cb = onProgress ? function (f) { onProgress("count", f); } : null;
    var t0 = nowMs();
    var value = null;
    var engine = "Lucy_Hedgehog/JS";
    if (x >= 1e7 && getWasmModule()) {
      value = primeCountWasm(x, cb);
      if (value !== null) engine = "compiled C/WebAssembly";
    }
    if (value === null) value = primeCount(x, cb);
    if (onProgress) onProgress("done", 1);
    return { x: x, value: value, engine: engine, ms: round1(nowMs() - t0) };
  }

  // ------------------------------------------------------------------
  // n-th prime via count + walk
  // ------------------------------------------------------------------

  var WINDOW_MAX = 1 << 23; // cap on integers per sieve window during the walk

  // Size the first walk window from the known deficit: |n − π(x0)| primes
  // at an average gap of ln x0, oversized 8×.  Later windows double, so a
  // pessimal estimate costs only a few extra rounds.
  function walkWindow(gapPrimes, x0) {
    var w = Math.ceil((gapPrimes + 16) * Math.log(x0) * 8);
    if (w < 1 << 16) w = 1 << 16;
    if (w > WINDOW_MAX) w = WINDOW_MAX;
    return w;
  }

  function nthPrimeByCounting(n, onProgress) {
    var tEst = nowMs();
    var x0 = guessNthPrime(n);
    var tCount = nowMs();
    var cb = onProgress ? function (f) { onProgress("count", f); } : null;
    var c0 = null;
    var engineName = "Lucy_Hedgehog/JS";
    if (x0 >= 1e7 && getWasmModule()) {
      c0 = primeCountWasm(x0, cb);
      if (c0 !== null) engineName = "compiled C/WebAssembly";
    }
    if (c0 === null) c0 = primeCount(x0, cb);
    return nthPrimeFromCount(n, x0, c0, {
      engine: engineName,
      msEstimate: round1(tCount - tEst),
      msCount: round1(nowMs() - tCount),
      onProgress: onProgress
    });
  }

  // Base primes for the walk cover √(any point it can visit): the n-th prime
  // is ≤ upperBoundForNthPrime(n) by Rosser's theorem.  Cached (a superset
  // is fine) so prepareWalk() can sieve them while the count runs elsewhere.
  var walkCache = null;
  function walkBasePrimes(n) {
    var limit = isqrt(upperBoundForNthPrime(n)) + 1;
    if (walkCache && walkCache.limit >= limit) return walkCache.primes;
    var primes = primesUpTo(limit);
    walkCache = { limit: limit, primes: primes };
    return primes;
  }
  function prepareWalk(n) {
    walkBasePrimes(checkN(n));
  }

  // Finish an n-th prime query from an exact count c0 = π(x0): sieve-walk to
  // the answer and attach the local verification.  Public so an external
  // counter (the multi-core engine) can supply c0.
  function nthPrimeFromCount(n, x0, c0, opts) {
    opts = opts || {};
    n = checkN(n);
    var onProgress = opts.onProgress;
    var tWalk = nowMs();
    var basePrimes = walkBasePrimes(n);
    var value = -1;
    var walked = 0;
    var seg, i;

    var win = walkWindow(Math.abs(n - c0), x0);

    if (c0 >= n) {
      // p(n) ≤ x0: walk downward.  cAbove = π(hi) for the current window top.
      var hi = x0;
      var cAbove = c0;
      for (;;) {
        var lo = Math.max(3, hi - win + 1);
        seg = sieveWindow(lo, hi, basePrimes);
        var cseg = countZeros(seg.comp, seg.len);
        walked += hi - lo + 1;
        var below = cAbove - cseg; // π(lo − 1), counting the prime 2
        if (below < n) {
          // target is the (n − below)-th prime inside this window
          var want = n - below;
          for (i = 0; i < seg.len; i++) {
            if (seg.comp[i] === 0 && --want === 0) {
              value = seg.base + 2 * i;
              break;
            }
          }
          break;
        }
        cAbove = below;
        hi = lo - 1;
        win = Math.min(win * 2, WINDOW_MAX);
        if (hi < 2) throw new Error("internal error: walked below 2");
      }
    } else {
      // p(n) > x0: walk upward from x0 + 1.
      var cnt = c0;
      var wlo = x0 + 1;
      for (;;) {
        var whi = wlo + win - 1;
        seg = sieveWindow(Math.max(3, wlo), whi, basePrimes);
        var cw = countZeros(seg.comp, seg.len);
        walked += whi - wlo + 1;
        if (cnt + cw >= n) {
          var need = n - cnt;
          for (i = 0; i < seg.len; i++) {
            if (seg.comp[i] === 0 && --need === 0) {
              value = seg.base + 2 * i;
              break;
            }
          }
          break;
        }
        cnt += cw;
        wlo = whi + 1;
        win = Math.min(win * 2, WINDOW_MAX);
      }
    }
    if (value < 0) throw new Error("internal error: walk failed to locate prime");
    var tEnd = nowMs();
    if (onProgress) onProgress("done", 1);
    var msEst = opts.msEstimate || 0, msCnt = opts.msCount || 0;
    return withNeighbors({
      n: n,
      value: value,
      method: "R-inverse estimate + Lucy_Hedgehog exact count [" + (opts.engine || "external counter") + "] + segmented sieve walk",
      guess: x0,
      piAtGuess: c0,
      offBy: n - c0, // primes between guess and answer (sign = direction)
      walked: walked,
      msEstimate: round1(msEst),
      msCount: round1(msCnt),
      msWalk: round1(tEnd - tWalk),
      ms: round1(msEst + msCnt + (tEnd - tWalk))
    });
  }

  function nowMs() {
    if (typeof performance !== "undefined" && performance.now) return performance.now();
    return Date.now();
  }

  function round1(x) {
    return Math.round(x * 10) / 10;
  }

  // ------------------------------------------------------------------
  // public API
  // ------------------------------------------------------------------

  function checkN(n) {
    n = typeof n === "string" ? Number(n) : n;
    if (typeof n !== "number" || !Number.isFinite(n) || !Number.isInteger(n)) {
      throw new RangeError("n must be an integer");
    }
    if (n < 1) throw new RangeError("n must be ≥ 1 (the 1st prime is 2)");
    if (n > MAX_N) throw new RangeError("n must be ≤ 2×10^14 (answers stay below 2^53)");
    return n;
  }

  // Promise-returning variant: opts.countAsync(x, onProgress) → Promise<π(x)>
  // (e.g. the multi-core engine in parallel.js) supplies the count; the
  // estimate and the walk are the same as nthPrime.
  function nthPrimeAsync(n, opts) {
    opts = opts || {};
    var n0;
    try { n0 = checkN(n); } catch (e) { return Promise.reject(e); }
    var minX = opts.parallelMinX || 1e12; // measured crossover vs the single-thread wasm core
    if (!opts.countAsync || n0 <= SIEVE_PATH_MAX) return Promise.resolve(nthPrime(n0, opts));
    var tEst = nowMs();
    var x0 = guessNthPrime(n0);
    if (x0 < minX) return Promise.resolve(nthPrime(n0, opts));
    var tCount = nowMs();
    var prog = opts.onProgress;
    var counting = opts.countAsync(x0, prog ? function (f) { prog("count", f); } : null);
    try { prepareWalk(n0); } catch (e) { /* the walk will sieve them itself */ }
    return counting.then(function (c0) {
      return nthPrimeFromCount(n0, x0, c0, {
        engine: opts.engineLabel || "external counter",
        msEstimate: round1(tCount - tEst),
        msCount: round1(nowMs() - tCount),
        onProgress: prog
      });
    });
  }

  function nthPrime(n, opts) {
    n = checkN(n);
    var onProgress = opts && opts.onProgress;
    var t0 = nowMs();

    if (n <= FIRST_PRIMES.length) {
      return withNeighbors({ n: n, value: FIRST_PRIMES[n - 1], method: "lookup table", ms: round1(nowMs() - t0) });
    }
    if (n <= SIEVE_PATH_MAX) {
      var v = nthPrimeBySieve(n);
      return withNeighbors({ n: n, value: v, method: "sieve of Eratosthenes (direct)", ms: round1(nowMs() - t0) });
    }
    return nthPrimeByCounting(n, onProgress); // neighbours attached inside
  }

  return {
    nthPrime: nthPrime,
    nthPrimeAsync: nthPrimeAsync,
    nthPrimeFromCount: nthPrimeFromCount,
    prepareWalk: prepareWalk,
    guessNthPrime: guessNthPrime,
    isPrime: isPrime,
    countPrimes: countPrimes,
    primeNeighbors: primeNeighbors,
    MAX_PI_X: MAX_PI_X,
    primeCount: primeCount,
    primeCountLMO: primeCountLMO,
    primeCountWasm: primeCountWasm,
    primeCountAuto: primeCountAuto,
    wasmAvailable: function () { return !!getWasmModule(); },
    icbrt: icbrt,
    riemannR: riemannR,
    nthPrimeEstimate: nthPrimeEstimate,
    inverseRiemannR: inverseRiemannR,
    primesUpTo: primesUpTo,
    upperBoundForNthPrime: upperBoundForNthPrime,
    lowerBoundForNthPrime: lowerBoundForNthPrime,
    isqrt: isqrt,
    MAX_N: MAX_N,
    SIEVE_PATH_MAX: SIEVE_PATH_MAX,
    version: "1.0.0"
  };
});
