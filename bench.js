/* Benchmark:  node bench.js [--big] [--huge]
 * --big adds n = 10^11, 10^12 (seconds); --huge adds 10^13, 10^14 (minutes). */
"use strict";

var NP = require("./nthprime.js");

var ns = [1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10];
if (process.argv.indexOf("--big") !== -1 || process.argv.indexOf("--huge") !== -1) ns.push(1e11, 1e12);
if (process.argv.indexOf("--huge") !== -1) ns.push(1e13, 1e14);

console.log("warming up JIT…");
NP.nthPrime(2e6);

console.log("\n        n |               p(n) |  time  | scaling");
console.log("----------+--------------------+--------+--------");
var prev = null;
ns.forEach(function (n) {
  var t = Date.now();
  var r = NP.nthPrime(n);
  var ms = Date.now() - t;
  var scale = "";
  if (prev && ms > 30 && prev.ms > 30) {
    // measured exponent e in time ∝ n^e between consecutive tiers
    var e = Math.log(ms / prev.ms) / Math.log(n / prev.n);
    scale = "n^" + e.toFixed(2);
  }
  console.log(
    pad("10^" + Math.round(Math.log10(n)), 9) + " | " +
    pad(String(r.value), 18) + " | " +
    pad(ms >= 1000 ? (ms / 1000).toFixed(1) + " s" : ms + " ms", 6) + " | " +
    pad(scale, 6)
  );
  prev = { n: n, ms: ms };
});
console.log("\n(exponent < 1.00 = sublinear: doubling n less than doubles the time)");

// ---- engine comparison on the largest tier that fits a coffee break:
// the compiled Deléglise–Rivat engine (single thread, the default) against
// the compiled Lucy_Hedgehog tables, single-thread and on every core
var PAR = null;
try { PAR = require("./parallel.js"); } catch (e) { PAR = null; }
var xCmp = ns.indexOf(1e12) !== -1 ? 3e13 : 3e11;
console.log("\nexact π(" + xCmp.toExponential(0).replace("+", "") + ") by engine:");
var tDR = Date.now(); var vDR = NP.primeCountDR(xCmp); tDR = Date.now() - tDR;
console.log("  Deléglise–Rivat (compiled, 1 thread) : " + (tDR / 1000).toFixed(2) + " s");
if (NP.wasmAvailable()) {
  var tL = Date.now(); var vL = NP.primeCountWasm(xCmp); tL = Date.now() - tL;
  console.log("  Lucy_Hedgehog   (compiled, 1 thread) : " + (tL / 1000).toFixed(2) + " s  " + (vL === vDR ? "same answer ✓" : "*** MISMATCH ***"));
}
if (PAR && PAR.available() && PAR.threads() >= 2) {
  var K = PAR.threads();
  var t2 = Date.now();
  PAR.primeCountParallel(xCmp, { threads: K }).then(function (vP) {
    var tP = Date.now() - t2;
    console.log("  Lucy_Hedgehog   (" + K + " threads)          : " + (tP / 1000).toFixed(2) + " s  " + (vP === vDR ? "same answer ✓" : "*** MISMATCH ***"));
  });
} else {
  console.log("  (multi-core engine unavailable in this environment)");
}

function pad(s, w) {
  while (s.length < w) s = " " + s;
  return s;
}
