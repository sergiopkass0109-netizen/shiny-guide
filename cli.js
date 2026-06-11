#!/usr/bin/env node
/* Command-line interface:
 *   node cli.js <n> [--json]            n-th prime, 1 ≤ n ≤ 2×10^14
 *   node cli.js --pi <x> [--engine lucy|lmo|both]   exact π(x)
 * Numbers accept 1234567, 1,234,567, 1_000_000, 1e9 and 10^12 forms. */
"use strict";

var NP = require("./nthprime.js");

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
var engine = "lucy";
var engIdx = args.indexOf("--engine");
if (engIdx !== -1) engine = args[engIdx + 1] || "lucy";
var positional = args.filter(function (a, i) {
  return a !== "--json" && a !== "--pi" && a !== "--engine" && i !== engIdx + 1;
});
var nArg = positional[0];

if (nArg === undefined) {
  console.error("usage: node cli.js <n> [--json]                   the n-th prime");
  console.error("       node cli.js --pi <x> [--engine lucy|lmo|both]   exact pi(x)");
  process.exit(2);
}

if (piMode) {
  var x = parseN(nArg);
  if (x === null || x < 0 || x > 9e15) {
    console.error("error: x must be a whole number up to 9×10^15 (got " + nArg + ")");
    process.exit(2);
  }
  var t0 = Date.now();
  if (engine === "both") {
    var va = NP.primeCount(x);
    var ta = Date.now() - t0;
    t0 = Date.now();
    var vb = NP.primeCountLMO(x);
    var tb = Date.now() - t0;
    console.log("pi(" + x.toLocaleString("en-US") + ") = " + va.toLocaleString("en-US") +
      "  [lucy " + ta + " ms]");
    console.log("pi(" + x.toLocaleString("en-US") + ") = " + vb.toLocaleString("en-US") +
      "  [lmo  " + tb + " ms]");
    console.log(va === vb ? "engines agree ✓" : "*** ENGINES DISAGREE — please report this x ***");
    process.exit(va === vb ? 0 : 1);
  }
  var val = engine === "lmo" ? NP.primeCountLMO(x) : NP.primeCount(x);
  console.log("pi(" + x.toLocaleString("en-US") + ") = " + val.toLocaleString("en-US") +
    "  [" + engine + " " + (Date.now() - t0) + " ms]");
  process.exit(0);
}

var n = parseN(nArg);
if (n === null || n < 1 || n > NP.MAX_N) {
  console.error("error: n must be a whole number between 1 and 2×10^14 (got " + nArg + ")");
  process.exit(2);
}

try {
  var res = NP.nthPrime(n);
  if (json) {
    console.log(JSON.stringify(res));
  } else {
    console.log("p(" + n.toLocaleString("en-US") + ") = " + res.value.toLocaleString("en-US"));
    console.log("  method : " + res.method);
    if (res.guess !== undefined) {
      console.log("  guess  : " + res.guess.toLocaleString("en-US") +
        "   (exact π there: " + res.piAtGuess.toLocaleString("en-US") +
        ", missed by " + Math.abs(res.offBy) + " primes)");
      console.log("  walk   : " + res.walked.toLocaleString("en-US") + " integers sieved");
    }
    console.log("  time   : " + res.ms + " ms" +
      (res.msCount !== undefined ? "  (count " + res.msCount + " ms, walk " + res.msWalk + " ms)" : ""));
  }
} catch (e) {
  console.error("error: " + ((e && e.message) || e));
  process.exit(1);
}
