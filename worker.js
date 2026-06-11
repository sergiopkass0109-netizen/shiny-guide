/* Web Worker: runs the nth-prime engine off the main thread. */
"use strict";
importScripts("nthprime.js");

self.onmessage = function (e) {
  var id = e.data.id;
  var n = e.data.n;
  try {
    var result = self.NthPrime.nthPrime(n, {
      onProgress: function (phase, frac) {
        self.postMessage({ id: id, progress: { phase: phase, frac: frac } });
      }
    });
    self.postMessage({ id: id, n: n, ok: true, result: result });
  } catch (err) {
    self.postMessage({ id: id, n: n, ok: false, error: String((err && err.message) || err) });
  }
};
