/*
 * Verification suite for nthprime.js
 *
 * Layers of defence:
 *   1. primesUpTo vs brute-force trial division (small range).
 *   2. nthPrime (table + sieve paths) vs an independently generated prime list.
 *   3. The counting path vs the sieve path on overlapping ranges — two
 *      independent algorithms must agree (incl. both sides of the cutoff).
 *   4. primeCount vs published values of π(10^k)   [Lehmer, OEIS A006880].
 *   5. nthPrime vs published values of the 10^k-th prime [OEIS A006988];
 *      p(10^9) = 22 801 763 489 and π(10^12) = 37 607 912 018 were
 *      re-confirmed against independent web sources for this project.
 *   6. Riemann R sanity: |R(10^k) − π(10^k)| must stay tiny (≪ √x).
 *   7. Estimate stays inside the rigorous Rosser/Dusart bracket.
 *   8. Input validation errors.
 *
 * Run:  node test/test.js          (~ a few seconds)
 *       node test/test.js --slow   (adds 10^10..10^12 anchors + big cross-engine checks; ~1 min)
 *       node test/test.js --huge   (adds p(10^13) and π(10^14), both engines; several minutes)
 */
"use strict";

var NP = require("../nthprime.js");

var SLOW = process.argv.indexOf("--slow") !== -1 || process.argv.indexOf("--huge") !== -1;
var HUGE = process.argv.indexOf("--huge") !== -1;
var failures = 0;
var checks = 0;

function check(cond, label) {
  checks++;
  if (!cond) {
    failures++;
    console.error("  FAIL  " + label);
  }
}

function eq(actual, expected, label) {
  check(actual === expected, label + "  (expected " + expected + ", got " + actual + ")");
}

function section(name) {
  console.log("\n== " + name);
}

function fmtMs(t0) {
  return " [" + (Date.now() - t0) + " ms]";
}

// ---------------------------------------------------------------- 1
section("primesUpTo vs trial division");
(function () {
  function isPrimeSlow(x) {
    if (x < 2) return false;
    for (var d = 2; d * d <= x; d++) if (x % d === 0) return false;
    return true;
  }
  var limit = 2000;
  var got = Array.from(NP.primesUpTo(limit));
  var want = [];
  for (var x = 2; x <= limit; x++) if (isPrimeSlow(x)) want.push(x);
  eq(got.length, want.length, "count of primes ≤ " + limit);
  for (var i = 0; i < want.length; i++) {
    if (got[i] !== want[i]) {
      check(false, "prime #" + (i + 1) + " mismatch: " + got[i] + " vs " + want[i]);
      break;
    }
  }
  check(true, "lists match");
  // edge limits
  eq(NP.primesUpTo(1).length, 0, "no primes ≤ 1");
  eq(NP.primesUpTo(2)[0], 2, "primes ≤ 2");
  eq(Array.from(NP.primesUpTo(3)).join(","), "2,3", "primes ≤ 3");
})();

// ---------------------------------------------------------------- 2
section("nthPrime table/sieve paths vs independent prime list");
(function () {
  var t0 = Date.now();
  var list = NP.primesUpTo(8000000); // 8e6 > p(500000) = 7368787
  for (var n = 1; n <= 2000; n++) {
    var v = NP.nthPrime(n).value;
    if (v !== list[n - 1]) {
      check(false, "nthPrime(" + n + ") = " + v + " ≠ " + list[n - 1]);
      return;
    }
  }
  check(true, "n = 1..2000 all match");
  // strided + spot checks across the whole sieve path, incl. the cutoff
  var spots = [49999, 100000, 123456, 250000, 499999, 500000];
  for (var n2 = 2097; n2 <= 500000; n2 += 4999) spots.push(n2);
  for (var i = 0; i < spots.length; i++) {
    var n3 = spots[i];
    eq(NP.nthPrime(n3).value, list[n3 - 1], "nthPrime(" + n3 + ")");
  }
  console.log("  ok" + fmtMs(t0));
})();

// ---------------------------------------------------------------- 3
section("counting path vs sieve path (independent algorithms agree)");
(function () {
  var t0 = Date.now();
  var list = NP.primesUpTo(NP.upperBoundForNthPrime(2000000));
  // both sides of the SIEVE_PATH_MAX cutoff
  var ns = [500001, 500002, 600000, 777777, 1000000, 1500000, 1999999];
  // deterministic pseudo-random extras
  var seed = 123456789;
  for (var k = 0; k < 25; k++) {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    ns.push(500001 + (seed % 1499999));
  }
  for (var i = 0; i < ns.length; i++) {
    var n = ns[i];
    var got = NP.nthPrime(n).value;
    var want = list[n - 1];
    if (got !== want) {
      check(false, "nthPrime(" + n + ") = " + got + " ≠ " + want);
    }
  }
  check(true, ns.length + " cross-checks agree");
  console.log("  ok" + fmtMs(t0));
})();

// ---------------------------------------------------------------- 4
section("primeCount vs published π(10^k)");
(function () {
  var PI = {
    1: 4,
    2: 25,
    3: 168,
    4: 1229,
    5: 9592,
    6: 78498,
    7: 664579,
    8: 5761455,
    9: 50847534,
    10: 455052511,
    11: 4118054813
  };
  if (SLOW) PI[12] = 37607912018; // re-confirmed via web sources, 2026-06
  for (var k = 1; k <= (SLOW ? 12 : 11); k++) {
    var t0 = Date.now();
    eq(NP.primeCount(Math.pow(10, k)), PI[k], "π(10^" + k + ")");
    if (k >= 10) console.log("  π(10^" + k + ")" + fmtMs(t0));
  }
  // small edge cases
  eq(NP.primeCount(1), 0, "π(1)");
  eq(NP.primeCount(2), 1, "π(2)");
  eq(NP.primeCount(3), 2, "π(3)");
  eq(NP.primeCount(4), 2, "π(4)");
  eq(NP.primeCount(10), 4, "π(10)");
})();

// ---------------------------------------------------------------- 5
section("nthPrime vs published 10^k-th primes (OEIS A006988)");
(function () {
  var P10 = {
    1: 29,
    2: 541,
    3: 7919,
    4: 104729,
    5: 1299709,
    6: 15485863,
    7: 179424673,
    8: 2038074743,
    9: 22801763489 // billionth prime, re-confirmed via web sources, 2026-06
  };
  if (SLOW) {
    P10[10] = 252097800623;
    P10[11] = 2760727302517;
    P10[12] = 29996224275833;
  }
  for (var k = 1; k <= (SLOW ? 12 : 9); k++) {
    var t0 = Date.now();
    var res = NP.nthPrime(Math.pow(10, k));
    eq(res.value, P10[k], "p(10^" + k + ")");
    if (k >= 8) {
      console.log(
        "  p(10^" + k + ") = " + res.value + fmtMs(t0) +
        (res.piAtGuess !== undefined
          ? "  guess off by " + res.offBy + " primes, walked " + res.walked
          : "")
      );
    }
  }
})();

// ---------------------------------------------------------------- 6
section("Riemann R stays close to π (validates Gram series + ζ)");
(function () {
  var PI = {
    4: 1229, 5: 9592, 6: 78498, 7: 664579, 8: 5761455,
    9: 50847534, 10: 455052511, 11: 4118054813
  };
  for (var k = 4; k <= 11; k++) {
    var x = Math.pow(10, k);
    var diff = Math.abs(NP.riemannR(x) - PI[k]);
    // |π − R| is conjecturally O(√x / ln x); allow a very generous margin —
    // this still catches any real formula/precision bug instantly.
    var margin = 0.05 * Math.sqrt(x) + 10;
    check(diff < margin, "|R(10^" + k + ") − π| = " + diff.toFixed(1) + " < " + margin.toFixed(1));
  }
})();

// ---------------------------------------------------------------- 7
section("estimate respects the rigorous Rosser/Dusart bracket");
(function () {
  var ns = [1e6, 1234567, 1e7, 5e7, 1e8, 1e9, 1e10, 1e11, 1e12];
  for (var i = 0; i < ns.length; i++) {
    var n = ns[i];
    var lo = NP.lowerBoundForNthPrime(n);
    var hi = NP.upperBoundForNthPrime(n);
    var x = Math.round(NP.inverseRiemannR(n));
    check(lo <= x && x <= hi, "R⁻¹(" + n + ") = " + x + " inside [" + lo + ", " + hi + "]");
  }
})();

// ---------------------------------------------------------------- 8
section("input validation");
(function () {
  function throws(fn, label) {
    try {
      fn();
      check(false, label + " did not throw");
    } catch (e) {
      check(e instanceof RangeError, label + " throws RangeError");
    }
  }
  throws(function () { NP.nthPrime(0); }, "n = 0");
  throws(function () { NP.nthPrime(-5); }, "n = -5");
  throws(function () { NP.nthPrime(1.5); }, "n = 1.5");
  throws(function () { NP.nthPrime(NaN); }, "n = NaN");
  throws(function () { NP.nthPrime(Infinity); }, "n = Infinity");
  throws(function () { NP.nthPrime(2e14 + 1); }, "n = 2×10^14 + 1");
  eq(NP.nthPrime(1).value, 2, "n = 1");
  eq(NP.nthPrime(2).value, 3, "n = 2");
  eq(NP.nthPrime(12).value, 37, "n = 12 (last table entry)");
  eq(NP.nthPrime(13).value, 41, "n = 13 (first sieve entry)");
  eq(NP.nthPrime("1000000").value, 15485863, "string input accepted");
})();

// ---------------------------------------------------------------- 9
section("LMO engine agrees with Lucy_Hedgehog engine (independent algorithms)");
(function () {
  var t0 = Date.now();
  // structured + random x, several α values — two fundamentally different
  // counting algorithms must produce identical results everywhere
  var xs = [89, 97, 100, 1009, 65536, 1e6, 1e6 + 3, 9999991, 1e8];
  var seed = 987654321;
  for (var i = 0; i < 20; i++) {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    xs.push(1000 + (seed % 999000000));
  }
  for (i = 0; i < xs.length; i++) {
    var x = xs[i];
    var a = NP.primeCount(x);
    var b = NP.primeCountLMO(x);
    if (a !== b) check(false, "engines disagree at x=" + x + ": lucy=" + a + " lmo=" + b);
    if (i % 4 === 0) {
      var c = NP.primeCountLMO(x, null, { alpha: 2 + (i % 5) });
      if (a !== c) check(false, "lmo alpha variant disagrees at x=" + x);
    }
  }
  check(true, xs.length + " cross-engine checks agree");
  eq(NP.primeCountLMO(1e10), 455052511, "LMO π(10^10)");
  if (NP.wasmAvailable()) {
    for (i = 0; i < xs.length; i += 2) {
      var xw = xs[i];
      var aw = NP.primeCount(xw);
      var bw = NP.primeCountWasm(xw);
      if (aw !== bw) check(false, "wasm disagrees at x=" + xw + ": js=" + aw + " wasm=" + bw);
    }
    eq(NP.primeCountWasm(1e10), 455052511, "WASM π(10^10)");
    check(true, "compiled engine agrees (triple-engine verification)");
  } else {
    console.log("  (WebAssembly engine unavailable here — skipped)");
  }
  console.log("  ok" + fmtMs(t0));
  if (SLOW) {
    var t1 = Date.now();
    eq(NP.primeCountLMO(1e12), 37607912018, "LMO π(10^12)");
    eq(NP.primeCountLMO(1e13), 346065536839, "LMO π(10^13)");
    eq(NP.primeCount(1e13), 346065536839, "Lucy π(10^13)");
    if (NP.wasmAvailable()) eq(NP.primeCountWasm(1e13), 346065536839, "WASM π(10^13)");
    console.log("  big anchors ok" + fmtMs(t1));
  }
  if (HUGE) {
    var t2 = Date.now();
    eq(NP.nthPrime(1e13).value, 323780508946331, "p(10^13)  [OEIS A006988]");
    console.log("  p(10^13) ok" + fmtMs(t2));
    t2 = Date.now();
    eq(NP.primeCount(1e14), 3204941750802, "Lucy π(10^14)");
    eq(NP.primeCountLMO(1e14), 3204941750802, "LMO π(10^14)");
    console.log("  π(10^14) both engines ok" + fmtMs(t2));
  }
})();

// ---------------------------------------------------------------- 10
section("deterministic polynomial-time primality test (isPrime)");
(function () {
  // exhaustive agreement with trial division on [0, 2000]
  function slow(x) {
    if (x < 2) return false;
    for (var d = 2; d * d <= x; d++) if (x % d === 0) return false;
    return true;
  }
  for (var x = 0; x <= 2000; x++) {
    if (NP.isPrime(x).prime !== slow(x)) { check(false, "isPrime(" + x + ") wrong"); return; }
  }
  check(true, "matches trial division on 0..2000");
  // Carmichael numbers fool Fermat tests — must not fool Miller–Rabin
  var carmichael = [561, 1105, 1729, 2465, 41041, 825265];
  for (var i = 0; i < carmichael.length; i++) {
    check(NP.isPrime(carmichael[i]).prime === false, "Carmichael " + carmichael[i] + " is composite");
  }
  // primes this project itself computed and anchored against OEIS
  check(NP.isPrime(22801763489).prime === true, "p(10^9) is prime");
  check(NP.isPrime(29996224275833).prime === true, "p(10^12) is prime");
  check(NP.isPrime("3475385758524527").prime === true, "p(10^14) is prime (string input)");
  // Mersenne prime 2^61 − 1, far beyond the calculator range
  check(NP.isPrime("2305843009213693951").prime === true, "2^61−1 is prime");
  check(NP.isPrime("2305843009213693953").prime === false, "2^61+1 is composite");
  check(NP.isPrime(15485865).factor === "3", "factor reporting (smallest factor)");
  check(NP.isPrime(1).prime === false, "1 is not prime");
  function throwsRange(fn, label) {
    try { fn(); check(false, label + " did not throw"); }
    catch (e) { check(e instanceof RangeError, label + " throws RangeError"); }
  }
  throwsRange(function () { NP.isPrime("3400000000000000000000000"); }, "beyond 3.3e24");
})();

// ---------------------------------------------------------------- 11
section("countPrimes + primeNeighbors (the page's pi(x) mode and local verification)");
(function () {
  var list = NP.primesUpTo(2000000);
  var xs = [0, 1, 2, 3, 10, 97, 1e6, 12345678, 1e9];
  for (var i = 0; i < xs.length; i++) {
    var r = NP.countPrimes(xs[i]);
    eq(r.value, NP.primeCount(xs[i]), "countPrimes(" + xs[i] + ")");
    check(typeof r.engine === "string" && r.ms >= 0, "countPrimes metadata for " + xs[i]);
  }
  // neighbours match an independently generated prime list exactly
  var idx = [0, 1, 2, 3, 100, 1000, 78497, 100000];
  for (i = 0; i < idx.length; i++) {
    var k = idx[i];
    var nb = NP.primeNeighbors(list[k]);
    eq(nb.prev, k === 0 ? null : list[k - 1], "prev of " + list[k]);
    eq(nb.next, list[k + 1], "next of " + list[k]);
  }
  // every nthPrime code path (table, sieve, counting) carries neighbours
  var ns = [1, 2, 12, 13, 1000, 140000];
  for (i = 0; i < ns.length; i++) {
    var res = NP.nthPrime(ns[i]);
    eq(res.prev, ns[i] === 1 ? null : list[ns[i] - 2], "nthPrime(" + ns[i] + ").prev");
    eq(res.next, list[ns[i]], "nthPrime(" + ns[i] + ").next");
  }
  var big = NP.nthPrime(600000); // counting path
  check(big.prev < big.value && big.value < big.next &&
        NP.isPrime(big.prev).prime && NP.isPrime(big.next).prime &&
        big.next - big.prev < 1000, "nthPrime(600000) neighbours are prime and near");
  function throwsRange(fn, label) {
    try { fn(); check(false, label + " did not throw"); }
    catch (e) { check(e instanceof RangeError, label + " throws RangeError"); }
  }
  NP.prepareWalk(700000);   // counting path (n > 500000) reuses the prepared base primes
  var expect = NP.primesUpTo(NP.upperBoundForNthPrime(700000))[599999];
  eq(NP.nthPrime(600000).value, expect, "walk uses the prepared base-prime cache");
  throwsRange(function () { NP.countPrimes(9e15 + 2); }, "countPrimes beyond 9e15");
  throwsRange(function () { NP.countPrimes(1.5); }, "countPrimes non-integer");
})();

// ---------------------------------------------------------------- 12
section("index.html is self-contained and in sync (run `npm run build` if not)");
(function () {
  var fs = require("fs");
  var path = require("path");
  var root = path.join(__dirname, "..");
  var html = fs.readFileSync(path.join(root, "index.html"), "utf8");
  var engine = fs.readFileSync(path.join(root, "nthprime.js"), "utf8").trim();
  var style = fs.readFileSync(path.join(root, "style.css"), "utf8").trim();
  var wasmjs = fs.readFileSync(path.join(root, "engine-wasm.js"), "utf8").trim();
  var parjs = fs.readFileSync(path.join(root, "parallel.js"), "utf8").trim();
  check(html.indexOf(parjs) !== -1, "index.html embeds current parallel.js");
  check(html.indexOf(engine) !== -1, "index.html embeds current nthprime.js");
  check(html.indexOf(style) !== -1, "index.html embeds current style.css");
  check(html.indexOf(wasmjs) !== -1, "index.html embeds current engine-wasm.js");
  check(html.indexOf("/*BUILD:") === -1, "no unreplaced build markers");
  check(html.indexOf("<script src=") === -1 && html.indexOf("<link rel=\"stylesheet\"") === -1,
    "no external file references");
})();

// ---------------------------------------------------------------- 13
section("multi-core engine (SharedArrayBuffer + worker_threads) agrees exactly");
var parallelDone = (function () {
  var PAR;
  try { PAR = require("../parallel.js"); } catch (e) { PAR = null; }
  if (!PAR || !PAR.available()) { console.log("  (multi-core engine unavailable here — skipped)"); return Promise.resolve(); }
  var t0 = Date.now();
  // [x, threads]  — 0 = auto; small x exercises the coordinator-only path,
  // larger x the banded parallel updates; a prime x and an odd thread count
  // probe the band boundaries
  var cases = [[1e6, 2], [2e9, 3], [3000000007, 5], [1e10, 0], [12345678901, 4]];
  var chain = Promise.resolve();
  cases.forEach(function (c) {
    chain = chain.then(function () {
      return PAR.primeCountParallel(c[0], { threads: c[1] || undefined }).then(function (v) {
        eq(v, NP.primeCount(c[0]), "multi-core π(" + c[0] + ") threads=" + (c[1] || "auto"));
      });
    });
  });
  chain = chain.then(function () {
    return NP.nthPrimeAsync(1e9, {
      engineLabel: "multi-core (test)",
      parallelMinX: 1e9,
      countAsync: function (x, cb) { return PAR.primeCountParallel(x, { onProgress: cb }); }
    }).then(function (r) {
      eq(r.value, 22801763489, "nthPrimeAsync(10^9) via multi-core = billionth prime");
      check(r.prev < r.value && r.value < r.next, "…carries verified neighbours");
      check(/multi-core/.test(r.method), "…labels the engine");
    });
  });
  chain = chain.then(function () {
    var job = PAR.primeCountParallel(1e11, { threads: 2 });
    job.cancel();
    return job.then(function () { check(false, "cancelled job resolved"); },
                    function (e) { check(/cancelled/.test(e.message), "cancel() rejects the job promptly"); });
  });
  return chain.then(function () { console.log("  ok (" + PAR.threads() + " threads here)" + fmtMs(t0)); },
                    function (err) { check(false, "multi-core engine error: " + (err && err.message)); });
})();

parallelDone.then(function () {
  console.log("\n" + checks + " checks, " + failures + " failure(s)" + (SLOW ? " [slow mode]" : ""));
  if (failures > 0) process.exit(1);
  console.log("ALL TESTS PASSED");
});
