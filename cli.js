#!/usr/bin/env node
/* Command-line interface:
 *   node cli.js <n> [--json] [--threads K]          n-th prime, 1 ≤ n ≤ 2×10^14 (multi-core when K ≠ 1)
 *   node cli.js --pi <x> [--engine lucy|lmo|wasm|parallel|all]   exact π(x)
 * Numbers accept 1234567, 1,234,567, 1_000_000, 1e9 and 10^12 forms. */
"use strict";

var NP = require("./nthprime.js");
var PAR = null;
try { PAR = require("./parallel.js"); } catch (e) { PAR = null; }

function parseN(raw) {
  var s = String(raw).trim().toLowerCase().replace(/[,_\s]/g, "");
  if (!s) return null;
  var v;
  var m = /^10\^(\d{1,2})$/.exec(s);
  if (m) v = Math.pow(10, +m[1]);
  else if (/^\d+$/.test(s) || /^\d+(\.\d+)?e\+?\d{1,2}$/.test(s)) v = Number(s);
  else return null;
  if (!isFinite(v) || Math.floor(v) !== v) return null;
  return v;
}

var args = process.argv.slice(2);
var json = args.indexOf("--json") !== -1;
var piMode = args.indexOf("--pi") !== -1;
var engIdx = args.indexOf("--engine");
var engine = engIdx !== -1 ? (args[engIdx + 1] || "lucy") : "lucy";
var thrIdx = args.indexOf("--threads");
var threads = thrIdx !== -1 ? Number(args[thrIdx + 1]) || 0 : 0;
var positional = args.filter(function (a, i) {
  return a !== "--json" && a !== "--pi" && a !== "--engine" && a !== "--threads" &&
    (engIdx === -1 || i !== engIdx + 1) && (thrIdx === -1 || i !== thrIdx + 1);
});
var nArg = positional[0];

if (nArg === undefined) {
  console.error("usage: node cli.js <n> [--json] [--threads K]                        the n-th prime");
  console.error("       node cli.js --pi <x> [--engine lucy|lmo|wasm|parallel|all]   exact pi(x)");
  process.exit(2);
}

var fmt = function (v) { return Number(v).toLocaleString("en-US"); };
var parallelOK = PAR && PAR.available() && threads !== 1;
var K = parallelOK ? (threads || PAR.threads()) : 1;

if (piMode) {
  var x = parseN(nArg);
  if (x === null || x < 0 || x > 9e15) {
    console.error("error: x must be a whole number up to 9×10^15 (got " + nArg + ")");
    process.exit(2);
  }
  var t0 = Date.now();
  if (engine === "both" || engine === "all") {
    var va = NP.primeCount(x);
    var ta = Date.now() - t0;
    t0 = Date.now();
    var vb = NP.primeCountLMO(x);
    var tb = Date.now() - t0;
    console.log("pi(" + fmt(x) + ") = " + fmt(va) + "  [lucy " + ta + " ms]");
    console.log("pi(" + fmt(x) + ") = " + fmt(vb) + "  [lmo  " + tb + " ms]");
    var allAgree = va === vb;
    var chain = Promise.resolve();
    if (engine === "all" && NP.wasmAvailable()) {
      t0 = Date.now();
      var vc = NP.primeCountWasm(x);
      console.log("pi(" + fmt(x) + ") = " + fmt(vc) + "  [wasm " + (Date.now() - t0) + " ms]");
      allAgree = allAgree && vc === va;
    }
    if (engine === "all" && parallelOK) {
      chain = chain.then(function () {
        var t1 = Date.now();
        return PAR.primeCountParallel(x, { threads: K }).then(function (vd) {
          console.log("pi(" + fmt(x) + ") = " + fmt(vd) + "  [parallel ×" + K + " " + (Date.now() - t1) + " ms]");
          allAgree = allAgree && vd === va;
        });
      });
    }
    chain.then(function () {
      console.log(allAgree ? "engines agree ✓" : "*** ENGINES DISAGREE — please report this x ***");
      process.exit(allAgree ? 0 : 1);
    });
  } else if (engine === "parallel") {
    if (!parallelOK) { console.error("error: multi-core engine unavailable here"); process.exit(1); }
    PAR.primeCountParallel(x, { threads: K }).then(function (v) {
      console.log("pi(" + fmt(x) + ") = " + fmt(v) + "  [parallel ×" + K + " " + (Date.now() - t0) + " ms]");
    }, function (e) { console.error("error: " + e.message); process.exit(1); });
  } else {
    var val = engine === "lmo" ? NP.primeCountLMO(x)
            : engine === "wasm" ? NP.primeCountWasm(x)
            : NP.primeCount(x);
    console.log("pi(" + fmt(x) + ") = " + fmt(val) + "  [" + engine + " " + (Date.now() - t0) + " ms]");
  }
} else {
  var n = parseN(nArg);
  if (n === null || n < 1 || n > NP.MAX_N) {
    console.error("error: n must be a whole number between 1 and 2×10^14 (got " + nArg + ")");
    process.exit(2);
  }
  var report = function (res) {
    if (json) { console.log(JSON.stringify(res)); return; }
    console.log("p(" + fmt(n) + ") = " + fmt(res.value));
    console.log("  method : " + res.method);
    if (res.guess !== undefined) {
      console.log("  guess  : " + fmt(res.guess) + "   (exact π there: " + fmt(res.piAtGuess) +
        ", missed by " + Math.abs(res.offBy) + " primes)");
      console.log("  walk   : " + fmt(res.walked) + " integers sieved");
    }
    if (res.prev !== undefined) {
      console.log("  verify : prime ✓  previous " + fmt(res.prev) + " (gap " + (res.value - res.prev) +
        ")  next " + fmt(res.next) + " (gap " + (res.next - res.value) + ")");
    }
    console.log("  time   : " + res.ms + " ms" +
      (res.msCount !== undefined ? "  (count " + res.msCount + " ms, walk " + res.msWalk + " ms)" : ""));
  };
  var run = parallelOK
    ? NP.nthPrimeAsync(n, {
        engineLabel: "multi-core ×" + K + " threads",
        parallelMinX: Math.max(2e11, 1e12 * 4 / K),
        countAsync: function (x0, cb) { return PAR.primeCountParallel(x0, { threads: K, onProgress: cb }); }
      })
    : Promise.resolve(NP.nthPrime(n));
  run.then(report, function (e) { console.error("error: " + ((e && e.message) || e)); process.exit(1); });
}
