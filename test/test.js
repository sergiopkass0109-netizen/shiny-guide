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
 *       node test/test.js --slow   (adds n = 10^10, 10^11 and π(10^12); minutes)
 */
"use strict";

var NP = require("../nthprime.js");

var SLOW = process.argv.indexOf("--slow") !== -1;
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
  }
  for (var k = 1; k <= (SLOW ? 11 : 9); k++) {
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
  throws(function () { NP.nthPrime(1e12 + 1); }, "n = 10^12 + 1");
  eq(NP.nthPrime(1).value, 2, "n = 1");
  eq(NP.nthPrime(2).value, 3, "n = 2");
  eq(NP.nthPrime(12).value, 37, "n = 12 (last table entry)");
  eq(NP.nthPrime(13).value, 41, "n = 13 (first sieve entry)");
  eq(NP.nthPrime("1000000").value, 15485863, "string input accepted");
})();

// ----------------------------------------------------------------
console.log("\n" + checks + " checks, " + failures + " failure(s)" + (SLOW ? " [slow mode]" : ""));
if (failures > 0) process.exit(1);
console.log("ALL TESTS PASSED");
