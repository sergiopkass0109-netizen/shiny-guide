#!/usr/bin/env node
/* Command-line interface:  node cli.js <n> [--json]
 * <n> accepts 1234567, 1,234,567, 1_000_000, 1e9 and 10^12 forms. */
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
var nArg = args.filter(function (a) { return a !== "--json"; })[0];

if (nArg === undefined) {
  console.error("usage: node cli.js <n> [--json]      e.g.  node cli.js 1e9");
  process.exit(2);
}

var n = parseN(nArg);
if (n === null || n < 1 || n > NP.MAX_N) {
  console.error("error: n must be a whole number between 1 and 10^12 (got " + nArg + ")");
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
