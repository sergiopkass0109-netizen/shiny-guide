/* Benchmark:  node bench.js [--big] [--huge]
 * --big adds n = 10^11, 10^12 (seconds); --huge adds 10^13, 10^14 (minutes). */
"use strict";

var NP = require("./nthprime.js");

var ns = [1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10];
if (process.argv.indexOf("--big") !== -1 || process.argv.indexOf("--huge") !== -1) ns.push(1e11, 1e12);
if (process.argv.indexOf("--huge") !== -1) ns.push(1e13, 1e14);

console.log("warming up JIT…");
NP.nthPrime(2e6);

console.log("\n        n |               p(n) |  time");
console.log("----------+--------------------+--------");
ns.forEach(function (n) {
  var t = Date.now();
  var r = NP.nthPrime(n);
  var ms = Date.now() - t;
  console.log(
    pad("10^" + Math.round(Math.log10(n)), 9) + " | " +
    pad(String(r.value), 18) + " | " +
    pad(ms >= 1000 ? (ms / 1000).toFixed(1) + " s" : ms + " ms", 6)
  );
});

function pad(s, w) {
  while (s.length < w) s = " " + s;
  return s;
}
