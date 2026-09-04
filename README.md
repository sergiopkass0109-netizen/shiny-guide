# The *n*-th Prime

[![CI](https://github.com/sergiopkass0109-netizen/shiny-guide/actions/workflows/ci.yml/badge.svg)](https://github.com/sergiopkass0109-netizen/shiny-guide/actions/workflows/ci.yml)

Type a number **n** into the box, get back the **n-th prime number** — exactly,
for any **1 ≤ n ≤ 2×10¹⁴**. Counts with the **Deléglise–Rivat algorithm** — the
research-grade method behind every π(x) record, compiled from C to
WebAssembly, in O(√x) memory — runs as an **installable offline app**, and every
answer is checked by **four independent engines** plus an on-the-spot
primality proof, with a one-click cross-check by a second algorithm on every
core. Zero dependencies.

```
p(10⁶)    =             15,485,863       ~10 ms
p(10⁹)    =         22,801,763,489       ~25 ms     (the billionth prime)
p(10¹²)   =     29,996,224,275,833       0.7 s      (the trillionth prime)
p(10¹³)   =    323,780,508,946,331       2.9 s
p(10¹⁴)   =  3,475,385,758,524,527        13 s      (matches OEIS A006988)
p(2×10¹⁴) =  7,093,600,525,704,677        21 s      (the float64 frontier)

π(10¹⁵)   =     29,844,570,422,669       5.4 s      (exact prime count, one thread, a few MB)
```

*(Node 22, one thread, the compiled Deléglise–Rivat engine of version 2.8 —
run `npm run bench` yourself; in-browser times are similar on a desktop.
For scale: version 2.7 needed 95 s for p(10¹⁴) on one thread and 33 s on
four; 2.4 needed 5 minutes.)*

## Multi-core: every CPU core, provably exact

Since 2.8 the multi-core Lucy_Hedgehog engine is the **run-time cross-check**
(the *cross-check* button under a result recounts π at the guess on every
core with this second algorithm and compares) and the counting engine on
machines with 12+ threads; the single-thread Deléglise–Rivat engine is faster
than it on up to ~8 threads at every size we measured. How it stays exact
across threads is unchanged:

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

Measured on this 4-vCPU machine (version 2.8, all answers identical):

| exact π(x) | Deléglise–Rivat, 1 thread | Lucy, 1 thread | Lucy, 2 threads | Lucy, 4 threads |
|---|---|---|---|---|
| π(10¹²) | **0.15 s** | 0.38 s | 0.27 s | 0.21 s |
| π(10¹³) | **0.30 s** | 0.90 s | 0.70 s | 0.54 s |
| π(3×10¹³) | **0.57 s** | 1.9 s | 1.4 s | 0.95 s |
| π(10¹⁴) | **1.3 s** | 4.5 s | 2.9 s | 1.9 s |

The page uses the multi-core engine for counting only on machines with 12
or more threads (where it catches up with Deléglise–Rivat) and otherwise
keeps it for the cross-check. The
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
   Firefox on a desktop or laptop**. The count itself is single-threaded and
   needs nothing special; the multi-core *cross-check* needs cross-origin
   isolation, which the app's service worker provides there (the very first
   visit reloads once to turn it on; Safari doesn't honour the worker-added
   headers, so it has no cross-check). Phones are fine to ~10¹².
2. **Install it** (*install as app* in the header): it opens in its own
   window with the service worker already warm, and works offline.
3. **Plug in.** A laptop on battery may throttle the core doing the work.
4. **Know the sizes.** Up to 10¹² is under a second; 10¹³ about three
   seconds; 10¹⁴ about fifteen; the top of the range (2×10¹⁴) about twenty,
   all on one thread and a few megabytes of memory. Use *cancel* freely.
5. **The command line** runs the same engines: `node cli.js 2e14`, and
   `node cli.js --pi 1e15 --engine all` runs all five (Deléglise–Rivat, two
   Lucy builds, LMO, multi-core Lucy) and confirms they agree.
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
| **this project** | **2×10¹⁴** (π to 9×10¹⁵) | **your browser tab / a single HTML file** |
| [primecount](https://github.com/kimwalisch/primecount) (Kim Walisch) | ~10²⁴ | native C++, multithreaded, installed |

So: **200× beyond the reference web tool, exact, running client-side**, with
the same algorithm family primecount uses (Deléglise–Rivat) and every answer
reproducible by three independent engines. What this project does *not*
claim: beating primecount. Its Gourdon variant with AVX and OpenMP does
π(10¹⁵) in about a second on one core where this engine takes 5.4 s in
WebAssembly, and it scales to 10²⁹ with 128-bit arithmetic where this one
stops at 2⁵³. This is the research-grade mathematics brought to a place it
has never been deployed — a single self-contained HTML file you can save,
open offline, and audit — with an engine of a second, unrelated algorithm
one click away to confirm any answer.

## Scaling: sublinear, measured — and the polynomial-time question

**How does time grow as n grows?** Not exponentially — not even linearly.
The count step costs ~p(n)^⅔/log² ≈ (n ln n)^⅔/log², so a 10× larger n costs
about 4.5× more time, forever. Measured on this machine (version 2.8):

| step up | time ratio | linear would be | exponential would be |
|---|---|---|---|
| 10¹² → 10¹³ | 4.2× | 10× | astronomically worse |
| 10¹³ → 10¹⁴ | 4.5× | 10× | astronomically worse |
| 10¹⁴ → 2×10¹⁴ | 1.6× | 2× | — |

Run `npm run bench` and read the `scaling` column: it prints the measured
exponent e in *time ∝ nᵉ* between tiers — ≈ 0.65 now (0.78 with the Lucy
engine), i.e. **sublinear**. (In fact every named algorithm here is a polynomial-in-n
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

### 2. Count — exact π at the guess, on four engines

**Speed engine (2.8) — Deléglise–Rivat, compiled from C to WebAssembly.**
This is the algorithm family (Lagarias–Miller–Odlyzko → Deléglise–Rivat →
Gourdon) behind primecount and every π(x) record: O(x^⅔/log²x) time and
O(√x) memory, against O(x^¾) time and O(√x) *tables* for the Lucy engine
below. With y = α·∛x and a = π(y),

    π(x) = S₀ + S_special + a − 1 − P₂(x, a)

* **S₀** (ordinary leaves): Σ μ(n)·φ(x/n, 6) over n ≤ y with lpf(n) > 13,
  φ from the 30030-wheel table — O(y).
* **S_special**: every leaf (n, p_b) with n ≤ y < n·p_b, p_b < lpf(n)
  contributes −μ(n)·φ(⌊x/(n p_b)⌋, b−1). Deléglise–Rivat's split by the
  value v = ⌊x/(n p_b)⌋: **trivial** (v < p_b, φ = 1) counted in closed form
  from the π table; **easy** (v < p_b², φ = π(v) − b + 2) answered by one
  lookup in a fully sieved segment; **hard** (the rest) answered from the
  partial sieve state just before p_b is crossed off.
* **P₂**: Σ over y < p ≤ √x of π(x/p) − π(p) + 1, read from the same
  segments as x/p sweeps through them.

One segmented bit sieve over [1, x/y] (2²⁰ numbers per segment, in L2)
serves all three: the primes ≤ 13 come pre-removed by 64-bit wheel masks,
the primes 17…149 are crossed off word-wise with residue masks and a
popcount, larger primes multiple by multiple; per-block counters walked
monotonically answer each prime's hard leaves in amortised O(1); a prefix
popcount after the last prime answers the easy leaves and P₂ in O(1). The
whole thing is ~300 lines of C (`engine.c`), verified digit for digit against
the unchanged Lucy tables on hundreds of values and several α, and, at the
component level, against the JavaScript LMO engine (same φ(x, a) and P₂ for
the same y). α ≈ 8 is the measured optimum from 10¹¹ on (12 above 10¹⁴):
larger α shortens the sieve but multiplies the easy leaves.

Memory is the other half of the story: at x = 3×10¹³ the Lucy tables take
66 MB and at 9×10¹⁵ over 1 GB; Deléglise–Rivat uses a few MB — a prime
list up to √x, three tables up to y and one 128 KB segment — so the top of
the range no longer needs a desktop's worth of RAM.

**Second engine — Lucy_Hedgehog's algorithm compiled from C to WebAssembly**
(O(x^¾) time, O(√x) space): a dynamic programme over the ≤ 2√x distinct
values of ⌊x/k⌋ with the recurrence `S(v) ← S(v) − [S(⌊v/p⌋) − π(p−1)]`.
It is the multi-core engine (below) and the run-time cross-check; the same
C core (`engine.c`, zero libc) ships as ~15 KB wasm modules embedded base64
in the page. Version 2.5 made it **3.5× faster** than 2.4 at 10¹³ and above,
by profiling the compiled code loop by loop and changing only *how* the
same recurrence is evaluated:

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

**Reference engine — Lucy_Hedgehog in pure JavaScript**, the automatic
fallback wherever WebAssembly is unavailable, and the oracle every compiled
engine is tested against.

**Verification engine — Lagarias–Miller–Odlyzko (1985) in JavaScript**:
π(x) = φ(x,a) + a − 1 − P₂(x,a), φ split into ordinary and special leaves,
the special leaves answered by a segmented sieve with a Fenwick tree, P₂ in
one ascending sweep. It was the first complete LMO implementation to run
inside a web page, and it is the component-level reference for the compiled
Deléglise–Rivat engine: the suite requires the same φ(x, a) and P₂ from both
for the same y.

Four engines, two algorithm families, two languages: Deléglise–Rivat
computes your answer, and any of the others will tell you if it is wrong.
Run them all: `node cli.js --pi 1e12 --engine all`.

### 3. Walk — segmented sieve over the tiny gap

Exact π(guess) + the guess's tiny error ⇒ an odds-only segmented sieve of
Eratosthenes walks at most a few million integers to the exact n-th prime.
Milliseconds.

## How it is verified

`npm test` (~15 s) runs 13 independent layers, 280 checks (more in --slow / --huge):

1. `primesUpTo` vs brute-force trial division;
2. table/sieve paths vs an independently generated prime list;
3. counting path vs sieve path across the cutover and 25 pseudo-random n;
4. exact π(10^k), k = 1…11, vs published values (k ≤ 14 in `--huge`);
5. p(10^k), k = 1…9, vs [OEIS A006988][oeis] (k ≤ 13 in `--slow`/`--huge`) —
   p(10⁹) and π(10¹²) were re-confirmed against independent web sources;
6. |R(x) − π(x)| ≪ √x (guards the Gram series / ζ implementation);
7. estimates always inside the rigorous Rosser/Dusart bracket;
8. input-validation edge cases;
9. **four-engine agreement** (compiled Deléglise–Rivat vs compiled Lucy vs
   JS Lucy vs JS LMO) on structured + random x and several α values — the
   flagship check — plus 75 values chosen to cross the compiled Lucy engine's
   block-start threshold, segment boundaries and corrected prefix reads, and
   the Deléglise–Rivat components (φ(x,a), P₂, a) matched against the JS LMO
   engine for the same y;
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
| `nthprime.js` | engine selection, the JS engines, estimate + walk (inlined at build; Node-ready) |
| `build.js` | `npm run build` regenerates index.html |
| `engine.c` | the C core — Deléglise–Rivat counting and cache-blocked Lucy counting for wasm32 |
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
