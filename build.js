/* Generates index.html from index.template.html by inlining style.css and
 * nthprime.js, making the page one self-contained file.
 * Run:  node build.js   (or: npm run build)
 */
"use strict";

var fs = require("fs");
var path = require("path");

function read(f) {
  return fs.readFileSync(path.join(__dirname, f), "utf8");
}

var template = read("index.template.html");
var style = read("style.css").trim();
var engine = read("nthprime.js").trim();
var wasmjs = read("engine-wasm.js").trim();
var parjs = read("parallel.js").trim();

if (engine.indexOf("</script") !== -1 || style.indexOf("</style") !== -1 ||
    wasmjs.indexOf("</script") !== -1 || parjs.indexOf("</script") !== -1) {
  throw new Error("inlined source must not contain a closing tag literal");
}

// split/join keeps replacement text literal (no $-pattern interpretation)
var html = template
  .split("/*BUILD:STYLE*/").join(style)
  .split("/*BUILD:ENGINE*/").join(engine)
  .split("/*BUILD:WASM*/").join(wasmjs)
  .split("/*BUILD:PARALLEL*/").join(parjs);

if (html.indexOf("/*BUILD:") !== -1) {
  throw new Error("unreplaced BUILD marker left in output");
}

fs.writeFileSync(path.join(__dirname, "index.html"), html);
console.log("wrote index.html (" + (html.length / 1024).toFixed(1) + " KB, self-contained)");
