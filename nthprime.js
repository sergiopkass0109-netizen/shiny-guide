/*
 * nthprime.js — compute the n-th prime number, fast.
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

  var MAX_N = 1e12; // p(MAX_N) ≈ 3.0e13, comfortably below 2^53
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

  function primeCount(N, onProgress) {
    N = Math.floor(N);
    if (N < 2) return 0;
    if (N > 9e15) throw new RangeError("primeCount: x exceeds 2^53 safety bound");
    var r = isqrt(N);
    var small = new Float64Array(r + 1);
    var large = new Float64Array(r + 1); // large[i] = S(⌊N/i⌋)
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
      for (; i <= imax; i++) large[i] -= small[Math.floor(Np / i)] - sp1;
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
  // n-th prime via count + walk
  // ------------------------------------------------------------------

  var WINDOW_MAX = 1 << 21; // cap on integers per sieve window during the walk

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
    var c0 = primeCount(x0, onProgress ? function (f) { onProgress("count", f); } : null);
    var tWalk = nowMs();

    // Base primes cover √(any point the walk can visit): the n-th prime is
    // ≤ upperBoundForNthPrime(n) by Rosser's theorem.
    var basePrimes = primesUpTo(isqrt(upperBoundForNthPrime(n)) + 1);
    var value = -1;
    var walked = 0;
    var seg, i, num;

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
    return {
      n: n,
      value: value,
      method: "R-inverse estimate + Lucy_Hedgehog exact count + segmented sieve walk",
      guess: x0,
      piAtGuess: c0,
      offBy: n - c0, // primes between guess and answer (sign = direction)
      walked: walked,
      msEstimate: round1(tCount - tEst),
      msCount: round1(tWalk - tCount),
      msWalk: round1(tEnd - tWalk),
      ms: round1(tEnd - tEst)
    };
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

  function nthPrime(n, opts) {
    n = typeof n === "string" ? Number(n) : n;
    if (typeof n !== "number" || !Number.isFinite(n) || !Number.isInteger(n)) {
      throw new RangeError("n must be an integer");
    }
    if (n < 1) throw new RangeError("n must be ≥ 1 (the 1st prime is 2)");
    if (n > MAX_N) throw new RangeError("n must be ≤ 10^12 (answers stay below 2^53)");
    var onProgress = opts && opts.onProgress;
    var t0 = nowMs();

    if (n <= FIRST_PRIMES.length) {
      return { n: n, value: FIRST_PRIMES[n - 1], method: "lookup table", ms: round1(nowMs() - t0) };
    }
    if (n <= SIEVE_PATH_MAX) {
      var v = nthPrimeBySieve(n);
      return { n: n, value: v, method: "sieve of Eratosthenes (direct)", ms: round1(nowMs() - t0) };
    }
    return nthPrimeByCounting(n, onProgress);
  }

  return {
    nthPrime: nthPrime,
    primeCount: primeCount,
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
