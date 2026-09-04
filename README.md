# The *n*-th Prime

[![CI](https://github.com/sergiopkass0109-netizen/shiny-guide/actions/workflows/ci.yml/badge.svg)](https://github.com/sergiopkass0109-netizen/shiny-guide/actions/workflows/ci.yml)

Type a number **n** into the box, get back the **n-th prime number** — exactly,
for any **1 ≤ n ≤ 2×10¹⁴**. Counts on **every CPU core**, runs as an
**installable offline app**, and every answer is verified by **three independent
engines** plus an on-the-spot primality proof. Zero dependencies.

```
p(10⁶)    =             15,485,863       ~12 ms
p(10⁹)    =         22,801,763,489       ~50 ms     (the billionth prime)
p(10¹²)   =     29,996,224,275,833       ~4 s       (the trillionth prime)   2.0 s on 4 threads
p(10¹³)   =    323,780,508,946,331      ~21 s                                8.6 s on 4 threads
p(10¹⁴)   =  3,475,385,758,524,527       ~95 s      (matches OEIS A006988)   33 s on 4 threads
p(2×10¹⁴) =  7,093,600,525,704,677      ~150 s      (the float64 frontier)   53 s on 4 threads
```

*(Node 22, single thread unless noted, the WebAssembly core — run `npm run bench` yourself.
In-browser times are similar on a desktop machine. Version 2.5 made the
count 3.5× faster than 2.4 by cache-blocking the compiled engine; see Count.)*

## Multi-core: every CPU core, provably exact

Most of the work is the **block sweep** described under *Count* below, and
it has no data hazards at all: the threads simply claim 4096-entry chunks of
the range through one atomic counter until it is exhausted, so the
division-heavy low end of a range never lands on a single thread. The
remaining classical passes write `large[i]` while reading `large[i·p]` (a
higher index) and write `small[v]` while reading `small[⌊v/p⌋]` (a lower
index); processed in **geometric bands** — `[p^k, p^{k+1})` ascending for
`large`, `(r/p^{k+1}, r/p^k]` descending for `small` — every band reads only
indices outside itself that nothing has written yet in that pass, so each
band splits across threads with no locks in the hot loops and one barrier per
band or sweep. `parallel.js` runs that over one shared memory; the same
source is every thread, thread 0 coordinates, and the hot loops are the
**compiled C kernels** (`engine_par.c`, a thread-safe WebAssembly build of the
verified engine) with line-for-line JavaScript kernels as the fallback. In
Node it uses `worker_threads`; in the browser it needs cross-origin
isolation, which the app's service worker provides on GitHub Pages (one
automatic reload on the very first visit). The page says *multi-core ready ·
N threads* when it's on, and *wasm+simd kernels* next to the engine name when
the SIMD build is running.

Measured on this 4-vCPU machine (v2.5), answers identical to the single-thread core:

| | v2.4 single | v2.5 single | v2.5 2 threads | v2.5 4 threads |
|---|---|---|---|---|
| π(10¹²) | 1.38 s | 0.35–0.7 s | 0.48 s | 0.39 s |
| π(10¹³) | 6.0 s | 1.7–1.9 s | 1.5 s | 1.1–1.3 s |
| π(3×10¹³) | 12.8 s | 3.4–3.7 s | 2.5–2.8 s | 1.7–2.6 s |
| p(10¹²) end-to-end | ~9 s | 4.1 s | — | **2.0 s** |

(Ranges are repeat runs; this shared machine's timing noise is ±10–40%. The
2-thread column is what a typical laptop gets: π(3×10¹³), the count behind
p(10¹²), went from 8.9 s to 2.5 s there — a 3.5× cut.)

More cores, more speed. The page switches to multi-core at x ≈ 10¹² on 4
threads, earlier with more cores (below that, thread wake-ups cost more than
they save). While the cores count, the compute worker sieves the base primes
the final walk needs, so the walk takes ~0.1 s even at the top end. The
kernels ship in two builds — baseline and SIMD (auto-vectorised fills,
+10–20%) — and the loader picks the SIMD one only where
`WebAssembly.validate` accepts it, so old browsers keep working. The test
suite (2, 3, 4, 5 and auto threads, both kernel kinds, prime and composite
x, cancel, and repeated runs) requires digit-for-digit agreement with the
three single-thread engines — which is how a representation mismatch between
the C kernels and the JavaScript views (u64 vs f64 of the same bytes) was
caught before release.

**A barrier race, found and fixed in 2.6.** Stress-testing the multi-core
engine under CPU load (three busy processes on a 4-core box) produced a
wrong count about once in 40 runs — sometimes off by thousands, once by
10¹⁴. Logging every thread's executed generations showed one thread running
a task *twice* while another skipped it: the coordinator bumps the
generation counter and then notifies; a thread that had already seen the
new counter, finished its task and gone back to waiting was woken by that
late notify, re-read the same generation, ran the task again, and its second
decrement released the barrier before the slower thread had run it at all.
The fix is a re-check after waking (`if (gen === seen) continue`). 2.5 and
2.4 carried this bug; every check we ever ran in the browser was exact, but
a heavily loaded machine could have shown a wrong number. After the fix,
840 runs under the same load were all exact, and the suite now repeats
multi-core runs with both kernel kinds.

## Running it at full speed

1. **Use the hosted app** (the GitHub Pages link) in **Chrome, Edge or
   Firefox on a desktop or laptop**. Multi-core needs cross-origin isolation,
   which the app's service worker provides there; the header must read
   *multi-core ready · N threads*. The very first visit reloads once to turn
   it on. Safari doesn't honour the worker-added isolation headers, so it
   runs single-threaded; phones are fine to ~10¹¹.
2. **Install it** (*install as app* in the header): it opens in its own
   window with the service worker already warm, and works offline.
3. **Plug in, and close memory-hungry tabs.** The count is limited by memory
   bandwidth, not arithmetic, so a browser full of heavy tabs slows it; a
   laptop on battery may throttle cores.
4. **Know the sizes.** Up to 10¹² is a few seconds; 10¹³ is under ten seconds
   on 4 threads; 10¹⁴ is about half a minute on 4 threads (1½ minutes on one);
   anything ≥ 2×10¹³ needs ~1 GB of free RAM. Use *cancel* freely — it stops
   every thread at once.
5. **For the biggest jobs use the command line**, which has no browser memory
   ceiling: `node cli.js 2e14` uses every core (add `--threads K` to tune —
   the number of *physical* cores is usually best), and
   `node cli.js --pi 1e13 --engine all` runs all engines and confirms they agree.
6. **A saved copy of `index.html`** works anywhere, offline, single-threaded.

## It's an app, not a page

Served over HTTPS (the GitHub Pages link) the browser offers **Install** —
or press the *install as app* button in the header. It opens in its own
window, keeps its icon in your dock / start menu / home screen, and **works
fully offline** (the service worker caches the self-contained page). Saving
`index.html` still works too; that copy just runs single-threaded.

## Why 2×10¹⁴ is a hard wall — and why we drive right up to it

Every integer in this program is an IEEE-754 double. Doubles represent every
integer below 2⁵³ = 9,007,199,254,740,992 **exactly** — and the code carries a
proof (header of `nthprime.js`) that every floor-division it takes is exact in
that range. p(2×10¹⁴) ≈ 7.06×10¹⁵ is as close to that cliff as a power-of-ten-ish
cap safely gets. Going further means BigInt arithmetic, which is 10–50× slower —
so 2×10¹⁴ is the natural frontier of *fast* exact JavaScript. We compute to the
edge of it.

## How it compares to "the best ones out there" — honestly

| tool | reach | where it runs |
|---|---|---|
| typical online "nth prime" calculators | ~10⁶–10⁹ | server or browser |
| [The Nth Prime Page (t5k.org)](https://t5k.org/nthprime/) — the standard reference | 10¹² | **their server** |
| **this project** | **2×10¹⁴** | **your browser tab / a single HTML file** |
| [primecount](https://github.com/kimwalisch/primecount) (Kim Walisch) | ~10²⁴ | native C++, multithreaded, installed |

So: **200× beyond the reference web tool, exact, running client-side** — and
every answer is reproducible by two independent algorithms (below). What this
project does *not* claim: beating primecount. Hand-tuned, multithreaded C++
with decades of optimization is faster than any web page; this is the same
mathematics brought to a place it has never been deployed — a single
self-contained HTML file you can save, open offline, and audit.

## Scaling: sublinear, measured — and the polynomial-time question

**How does time grow as n grows?** Not exponentially — not even linearly.
The count step costs ~p(n)^¾ ≈ (n ln n)^¾, so a 10× larger n costs about
5.6× more time, forever. Measured on this machine:

| step up | time ratio | linear would be | exponential would be |
|---|---|---|---|
| 10¹² → 10¹³ | 5.1× | 10× | astronomically worse |
| 10¹³ → 10¹⁴ | 4.5× | 10× | astronomically worse |
| 10¹⁴ → 2×10¹⁴ | 1.6× | 2× | — |

Run `npm run bench` and read the `scaling` column: it prints the measured
exponent e in *time ∝ nᵉ* between tiers — consistently ≈ 0.78, i.e.
**sublinear**. (In fact every named algorithm here is a polynomial-in-n
algorithm; the field's open frontier is different — see below.)

**And "polynomial time" in the computer-science sense?** There the input
size is the number of *digits* d of n, and the honest state of mathematics
(checked against the 2026 literature) is:

* **Testing** whether a given number is prime in poly(d) time: **solved**
  (AKS 2002; deterministic Miller–Rabin below 3.3×10²⁴ — Sorenson & Webster
  2016). *Implemented in this project*: end any number with `?` in the box,
  or `NP.isPrime("4999999999999999999999")` — answers in O(d³), instantly,
  up to 3.3×10²⁴, a billion times beyond the calculator's own range.
* **Finding** the n-th prime (or even *any* d-digit prime, deterministically)
  in poly(d) time: **open problem** — the Polymath 4 project ("Finding
  primes", 2009–2012) attacked exactly this and it remains unsolved. The
  best known exact-π(x) theory (Lagarias–Odlyzko analytic, ~x^½) is still
  exponential in d. Anyone who solves it gets a very famous paper, and this
  README will be delighted to cite them.

## The mathematics (three stages + a twist)

### 1. Estimate — Riemann's R function, inverted

Cipolla's 1902 asymptotic expansion seeds Newton iteration on **Riemann's R
function**, computed via the [Gram series](https://mathworld.wolfram.com/GramSeries.html)
(ζ by Euler–Maclaurin):

```
R(x) = 1 + Σ_{k≥1} (ln x)^k / (k · k! · ζ(k+1))        x ← x − (R(x) − n)·ln x
```

Measured accuracy of the resulting guess for p(n):

| n | guess missed by |
|---|---|
| 10⁹ | 1,429 primes (6×10⁻⁵ %) |
| 10¹² | 35,884 primes (4×10⁻⁶ %) |
| 10¹⁴ | 49,262 primes (5×10⁻⁹ %) |

The guess is clamped into the **rigorous bracket**
n(ln n + ln ln n − 1) ≤ p(n) < n(ln n + ln ln n) (Dusart / Rosser), so a wild
estimate can never produce a wrong answer — only a slower one.

### 2. Count — exact π at the guess, on three engines

**Speed engine — Lucy_Hedgehog's algorithm compiled from C to WebAssembly**
(O(x^¾) time, O(√x) space): a dynamic programme over the ≤ 2√x distinct
values of ⌊x/k⌋ with the recurrence `S(v) ← S(v) − [S(⌊v/p⌋) − π(p−1)]`.
The C core (`engine.c`, zero libc) ships as ~5 KB wasm modules embedded
base64 in the page. Version 2.5 made it **3.5× faster** than 2.4 at 10¹³ and
above, by profiling the compiled code loop by loop and changing only *how*
the same recurrence is evaluated:

* **Runs, walked by quotient.** For i > √(x/p) consecutive i share one
  quotient q = ⌊x/(ip)⌋, and the run of q is exactly
  (⌊xp/(q+1)⌋, ⌊xp/q⌋]. Walking q downward costs one *independent*
  division per run; the previous form (`q = xp/i`, then `end = xp/q`, then
  `i = end+1`) was a serial chain of two divisions per run, ~30 cycles of
  pure latency — 52% of the whole count.
* **Blocks of primes, one cache-blocked sweep.** Primes are processed in
  blocks of up to 32 consecutive primes p₁ < … < pₖ ≤ 2p₁. For every
  i > I₀ = max(r/p₁, x/p₁³), each prime's update of `large[i]` reads only
  `small[q]` with q < p₁² — a region no small-pass of the block touches — so
  the k passes commute and are applied to `large[(I₀, imax]]` in one sweep
  through 32 KiB segments. Within a segment, each prime's runs go to a
  **difference array** (`D[start] += c; D[end+1] −= c`) and one prefix-sum
  pass materialises all k primes at once: the 3×10⁹ fill writes at x = 3×10¹³
  become 2×10⁸ writes into L1 cache. The prefix i ≤ I₀ is then done prime by
  prime in the classical order; its reads of `large[m]` with I₀ < m ≤ r add
  back the sweep's contributions of the later primes, so every value read is
  exactly S at the right stage.
* **Wheel start (2.7).** The tables begin at the state after the passes for
  2, 3, 5, 7, 11 and 13 — Lucy's invariant makes that state explicit:
  S(v) = #{m ≤ v coprime to 30030} − 1 + #{wheel primes ≤ v}, one table
  lookup per entry — instead of running those six full passes, which carried
  the largest heads of all (√(x/2) divisions for p = 2). Measured +6–9% in
  every mode.
* **Prefixes split around the sweep (2.6).** The prefix reads above I₀
  need one correction per later prime of the block; doing the first half of
  the block's prefixes *before* the sweep (correcting for the earlier primes
  instead) halves those ~4×10⁷ divisions at x = 3×10¹³ — measured +4–6% in
  every mode, single- and multi-thread.
* **What did *not* help, measured:** a reciprocal table instead of division
  (−26%), hand-written two-lane SIMD division (−17%, lane moves are
  expensive), signed float→int conversion (−44%), a split division/gather
  loop with `f64x2.div` (−14 to −21%). Auto-vectorised fills (+10–20%) and
  non-trapping float→int conversion (+13%) did, and stayed.
* **Index-major heads with per-prime reciprocals, rejected (2.7).** For one
  index the block's 32 quotients are ⌊y/p_j⌋ with the same y = ⌊x/i⌋, so one
  division plus 32 multiplications by 1/p_j (exact with a one-step
  correction, or with a cheap fraction test) could replace 32 divisions and
  32 read-modify-writes. Built, verified exact, and measured 0.62× and
  0.85×: in this runtime the loop is bound by instruction count, not by the
  divider, and the multiply-and-fix path issues more instructions than the
  division it removes. Block size (48, 64) and segment size (8192) changed
  nothing outside noise.
* **The primes above ∛x, three ways, all rejected (2.6).** For p³ > x the
  passes commute (every read is final), so they can be summed per index
  with no barriers at all. Three orders were built and verified exact:
  prime-major with per-thread accumulators (+13% single-thread here, but
  its reads scatter over the whole table and it stopped scaling across
  threads), pairs sorted by the product i·p in cache-sized windows (+2–6%
  single, −20% multi: the cursor bookkeeping cost more than the locality it
  bought), and groups of 32 primes per index (−3%). The block sweep already
  processes 32 primes per index segment — the same locality with no extra
  machinery — and stayed.

Every variant was accepted only after agreeing digit-for-digit with the
unchanged JavaScript engine on ~3400 inputs, then with the LMO engine in
the test suite. An older lesson still applies: this algorithm is bound by
64-bit division, and the exact trick in both languages is pipelined *double*
division with a proof of exactness below 2⁵³.

**Reference engine — the same algorithm in pure JavaScript**, used as the
automatic fallback wherever WebAssembly is unavailable, and as engine #2 in
verification.

**Verification engine #3 — Lagarias–Miller–Odlyzko (1985)**, the algorithm family
behind every prime-counting world record: π(x) = φ(x,a) + a − 1 − P₂(x,a),
with φ split into *ordinary leaves* (Möbius-weighted wheel counts) and
*special leaves* (answered by a segmented sieve with a Fenwick tree), and P₂
folded into one ascending segmented sweep. To our knowledge this is the
**first complete LMO implementation that runs inside a web page**.

Measured head-to-head, Lucy's constants beat LMO(α=1) ~2× at every size
below ~10¹⁷ — so Lucy computes your answer and LMO's job is *proof*: an
unrelated algorithm with an unrelated bug surface, agreeing digit-for-digit.
Run all three: `node cli.js --pi 1e12 --engine all`.

### 3. Walk — segmented sieve over the tiny gap

Exact π(guess) + the guess's tiny error ⇒ an odds-only segmented sieve of
Eratosthenes walks at most a few million integers to the exact n-th prime.
Milliseconds.

## How it is verified

`npm test` (~12 s) runs 12 independent layers, 260 checks (more in --slow / --huge):

1. `primesUpTo` vs brute-force trial division;
2. table/sieve paths vs an independently generated prime list;
3. counting path vs sieve path across the cutover and 25 pseudo-random n;
4. exact π(10^k), k = 1…11, vs published values (k ≤ 14 in `--huge`);
5. p(10^k), k = 1…9, vs [OEIS A006988][oeis] (k ≤ 13 in `--slow`/`--huge`) —
   p(10⁹) and π(10¹²) were re-confirmed against independent web sources;
6. |R(x) − π(x)| ≪ √x (guards the Gram series / ζ implementation);
7. estimates always inside the rigorous Rosser/Dusart bracket;
8. input-validation edge cases;
9. **triple-engine agreement** (compiled wasm vs JS Lucy vs LMO) on
   structured + random x and several α values — the flagship check — plus
   75 values chosen to cross the compiled engine's block-start threshold,
   segment boundaries and corrected prefix reads;
10. the deterministic polynomial-time primality test vs trial division,
    Carmichael numbers, Mersenne primes and the project's own verified primes;
11. `countPrimes` (the page's π(x) mode) and the neighbouring-prime
    verification attached to every answer, vs an independent prime list;
12. the self-contained `index.html` embeds the current engine byte-for-byte;
    and the multi-core engine, with both kernel kinds, on 2–5 and auto
    threads, including twelve repeated runs (the barrier regression).

Top-end results verified this way: π(10¹³) = 346,065,536,839 and
π(10¹⁴) = 3,204,941,750,802 (both engines, both correct), and
p(10¹³)/p(10¹⁴) match OEIS A006988 exactly.  p(2×10¹⁴) is beyond every
published table we know of — reproduce it yourself: `node cli.js 2e14`.

## Using it

**Web** — the live page, or open `index.html` (it is one self-contained file;
saving just it works):

```sh
python3 -m http.server 8000     # then http://localhost:8000
```

Three things the box understands:

| you type | you get |
|---|---|
| `1e12` (or `1,000,000,000,000`, `10^12`) | the n-th prime — with its neighbouring primes verified on the spot |
| `pi(1e12)` | π(x), the exact number of primes ≤ x, up to 9×10¹⁵ |
| `2305843009213693951?` | a deterministic primality verdict, up to 3.3×10²⁴ |

Every result gets a **share link** (e.g. `…/#n=1000000000000`) that reproduces
it on open, an elapsed-time readout, and a **cancel** button for the long ones.
Sizes ≥ 10¹² show a progress bar; n ≥ 2×10¹³ needs ~1 GB of memory — fine in
desktop Chrome/Firefox or the Node CLI, not for phones.

**CLI**

```sh
node cli.js 1e12                          # p(10^12) on every core (add --threads 1 for single)
node cli.js --pi 1e13 --engine all       # exact π(x) on all engines incl. multi-core, agreement check
```

**Library**

```js
const NP = require("./nthprime.js");
NP.nthPrime(1e12).value       // 29996224275833
NP.primeCount(1e12)           // 37607912018   (Lucy_Hedgehog)
NP.primeCountLMO(1e12)        // 37607912018   (Lagarias–Miller–Odlyzko)
```

**Tests / benchmarks**

```sh
npm test            # full fast suite (~10 s)
npm run test:slow   # adds 10^10..10^13 anchors, dual-engine (~2 min)
npm run test:huge   # adds p(10^13), π(10^14) both engines (~10 min)
npm run bench       # timing table; --big / --huge tiers
```

## Files

| file | purpose |
|---|---|
| `index.html` | the page — **generated**, fully self-contained |
| `index.template.html` | page source template |
| `style.css` | styling (inlined at build) |
| `nthprime.js` | both engines + estimate + walk (inlined at build; Node-ready) |
| `build.js` | `npm run build` regenerates index.html |
| `engine.c` | the C core — cache-blocked Lucy counting for wasm32 |
| `wasmbuild.js` | `npm run build:wasm` — clang → four .wasm builds → engine-wasm.js |
| `engine-wasm.js` | **generated** base64-embedded wasm loader with SIMD detection (committed so clang isn't required) |
| `engine.wasm` / `engine_simd.wasm` | the single-thread core, baseline and SIMD builds |
| `cli.js` | nth-prime and π(x) command line |
| `test/test.js` | the verification suite |
| `bench.js` | timing table |
| `parallel.js` | the multi-core engine (shared memory + workers; inlined into the page) |
| `engine_par.c` / `engine_par*.wasm` | thread-safe compiled kernels for the multi-core engine, baseline and SIMD builds (base64 in engine-wasm.js) |
| `sw.js` | service worker: offline cache + cross-origin isolation for multi-core |
| `manifest.webmanifest`, `icon-*.png`, `icon.js` | installable-app metadata; icons are generated by `node icon.js` |
| `.github/workflows/ci.yml` | runs the whole verification suite on every push (Node 20 + 22) |

## Research sources

* [primecount][primecount] — architecture (estimate → count → sieve) and the
  LMO/Deléglise–Rivat/Gourdon lineage
* Lagarias, Miller & Odlyzko, *Computing π(x): the Meissel–Lehmer method*,
  Math. Comp. 44 (1985) — the O(x^⅔) counting algorithm implemented here
* [Lucy_Hedgehog's method][lucy] (Project Euler #10, 2013; exposition by
  Jacob Elafandi and [griff's math blog](https://gbroxey.github.io/blog/2023/04/09/lucy-fenwick.html))
* [Gram series / Riemann prime counting — MathWorld][gram]
* Cipolla (1902) asymptotic for p(n); modern bounds:
  [Axler, *New estimates for the nth prime number* (2019)](https://cs.uwaterloo.ca/journals/JIS/VOL22/Axler/axler17.pdf),
  Dusart (1999), Rosser's theorem
* [OEIS A006988][oeis] (10^n-th primes) and
  [How many primes are there? (Prime Pages)](https://t5k.org/howmany.html)

[primecount]: https://github.com/kimwalisch/primecount
[lucy]: https://math.berkeley.edu/~elafandi/euler/p10/
[gram]: https://mathworld.wolfram.com/GramSeries.html
[oeis]: https://oeis.org/A006988
