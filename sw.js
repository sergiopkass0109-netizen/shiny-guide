/* Service worker with two jobs:
 *  1. offline: cache the app shell (the page is fully self-contained), network-first;
 *  2. multi-core: attach the cross-origin-isolation headers that static hosts like
 *     GitHub Pages cannot send, which is what unlocks SharedArrayBuffer for the
 *     multi-threaded engine.  No external resources are loaded, so the strict
 *     embedder policy costs nothing.
 */
var CACHE = "nthprime-v2.5.0";
var SHELL = ["./", "./index.html", "./manifest.webmanifest", "./icon-192.png", "./icon-512.png"];

self.addEventListener("install", function (e) {
  self.skipWaiting();
  e.waitUntil(caches.open(CACHE).then(function (c) { return c.addAll(SHELL); }).catch(function () {}));
});

self.addEventListener("activate", function (e) {
  e.waitUntil(
    caches.keys()
      .then(function (keys) { return Promise.all(keys.filter(function (k) { return k !== CACHE; }).map(function (k) { return caches.delete(k); })); })
      .then(function () { return self.clients.claim(); })
  );
});

self.addEventListener("fetch", function (e) {
  if (e.request.cache === "only-if-cached" && e.request.mode !== "same-origin") return;
  e.respondWith((async function () {
    var resp;
    try {
      resp = await fetch(e.request);
      if (resp && resp.ok && e.request.method === "GET" && new URL(e.request.url).origin === self.location.origin) {
        var c = await caches.open(CACHE);
        c.put(e.request, resp.clone());
      }
    } catch (err) {
      resp = (await caches.match(e.request)) || (await caches.match("./index.html"));
      if (!resp) throw err;
    }
    if (resp.status === 0) return resp; // opaque
    var h = new Headers(resp.headers);
    h.set("Cross-Origin-Embedder-Policy", "require-corp");
    h.set("Cross-Origin-Opener-Policy", "same-origin");
    return new Response(resp.body, { status: resp.status, statusText: resp.statusText, headers: h });
  })());
});
