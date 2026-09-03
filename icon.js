/* Generates icon-192.png and icon-512.png (a π glyph on the app's dark panel)
 * with nothing but Node's zlib — no image libraries.  Run: node icon.js */
"use strict";
var zlib = require("zlib");
var fs = require("fs");

var CRC = (function () {
  var t = new Uint32Array(256);
  for (var n = 0; n < 256; n++) {
    var c = n;
    for (var k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[n] = c >>> 0;
  }
  return t;
})();
function crc32(buf) {
  var c = 0xffffffff;
  for (var i = 0; i < buf.length; i++) c = CRC[(c ^ buf[i]) & 255] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}
function chunk(type, data) {
  var len = Buffer.alloc(4); len.writeUInt32BE(data.length);
  var td = Buffer.concat([Buffer.from(type, "ascii"), data]);
  var crc = Buffer.alloc(4); crc.writeUInt32BE(crc32(td));
  return Buffer.concat([len, td, crc]);
}
function png(size, pixel) {
  var stride = size * 4 + 1;
  var raw = Buffer.alloc(stride * size);
  for (var y = 0; y < size; y++) {
    raw[y * stride] = 0;
    for (var x = 0; x < size; x++) {
      var px = pixel(x / size, y / size);
      var o = y * stride + 1 + x * 4;
      raw[o] = px[0]; raw[o + 1] = px[1]; raw[o + 2] = px[2]; raw[o + 3] = 255;
    }
  }
  var ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(size, 0); ihdr.writeUInt32BE(size, 4);
  ihdr[8] = 8; ihdr[9] = 6; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
  return Buffer.concat([
    Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
    chunk("IHDR", ihdr), chunk("IDAT", zlib.deflateSync(raw, { level: 9 })), chunk("IEND", Buffer.alloc(0))
  ]);
}
// π made of bars: top bar + two stems; a faint diagonal glow behind it
function pixel(u, v) {
  var bg = [0x11, 0x16, 0x1f];
  var glow = Math.max(0, 1 - Math.hypot(u - 0.7, v - 0.25) / 0.9) * 0.25;
  var r = bg[0] + glow * 0x10, g = bg[1] + glow * 0x20, b = bg[2] + glow * 0x30;
  var inBar = v >= 0.30 && v <= 0.385 && u >= 0.21 && u <= 0.79;
  var inStem = v > 0.385 && v <= 0.73 && ((u >= 0.30 && u <= 0.385) || (u >= 0.615 && u <= 0.70));
  var foot = v > 0.73 && v <= 0.76 && ((u >= 0.27 && u <= 0.415) || (u >= 0.585 && u <= 0.73));
  if (inBar || inStem || foot) return [0x7e, 0xe7, 0x87];
  return [r | 0, g | 0, b | 0];
}
[192, 512].forEach(function (s) { fs.writeFileSync("icon-" + s + ".png", png(s, pixel)); });
console.log("wrote icon-192.png, icon-512.png");
