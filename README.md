# The *n*-th Prime

[![CI](https://github.com/sergiopkass0109-netizen/shiny-guide/actions/workflows/ci.yml/badge.svg)](https://github.com/sergiopkass0109-netizen/shiny-guide/actions/workflows/ci.yml)

Type a number **n** into the box, get back the **n-th prime number** — exactly,
for any **1 ≤ n ≤ 2×10¹⁴**, in a self-contained web page with zero dependencies —
counting on a **compiled C/WebAssembly core** at near-native speed, verified by
**three independent engines**.

```
p(10⁶)    =             15,485,863       ~10 ms
p(10⁹)    =         22,801,763,489       ~90 ms     (the billionth prime)
p(10¹²)   =     29,996,224,275,833       ~9 s       (the trillionth prime)
p(10¹³)   =    323,780,508,946,331      ~55 s
p(10¹⁴)   =  3,475,385,758,524,527       ~5 min     (matches OEIS A006988)
p(2×10¹⁴) =  7,093,600,525,704,677       ~9 min     (the float64 frontier)
```

*(Node 22, single thread, JavaScript with the WebAssembly core — run `npm run bench` yourself.
In-browser times are similar on a desktop machine.)*

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
| 10¹² → 10¹³ | 6.1× | 10× | astronomically worse |
| 10¹³ → 10¹⁴ | 5.5× | 10× | astronomically worse |
| 10¹⁴ → 2×10¹⁴ | 1.8× | 2× | — |

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
The C core (`engine.c`, ~60 lines, zero libc) ships as a 1.8 KB wasm module
embedded base64 in the page; steady-state it runs ~1.5× faster than the JS
engine. A performance lesson we measured the hard way: a naive C port was
*slower* than JavaScript (0.7×), because this algorithm is bound by 64-bit
integer division (~25–40 cycles each). The fix — in both languages — is
pipelined *double* division with a proof of exactness below 2⁵³. V8's JIT
had been applying that trick all along; C had to be taught it.

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

`npm test` (~10 s) runs 12 independent layers, 239 checks (more in --slow / --huge):

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
   structured + random x and several α values — the flagship check;
10. the deterministic polynomial-time primality test vs trial division,
    Carmichael numbers, Mersenne primes and the project's own verified primes;
11. `countPrimes` (the page's π(x) mode) and the neighbouring-prime
    verification attached to every answer, vs an independent prime list;
12. the self-contained `index.html` embeds the current engine byte-for-byte.

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
node cli.js 1e12                          # p(10^12) with diagnostics
node cli.js --pi 1e13 --engine both      # exact π(x), both engines, agreement check
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
| `engine.c` | the C core — Lucy counting for wasm32 |
| `wasmbuild.js` | `npm run build:wasm` — clang → engine.wasm → engine-wasm.js |
| `engine-wasm.js` | **generated** base64-embedded wasm loader (committed so clang isn't required) |
| `cli.js` | nth-prime and π(x) command line |
| `test/test.js` | the verification suite |
| `bench.js` | timing table |
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
