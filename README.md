# The *n*-th Prime

Type a number **n** into the box, get back the **n-th prime number** — exactly,
for any **1 ≤ n ≤ 2×10¹⁴**, in a self-contained web page with zero dependencies.

```
p(10⁶)    =             15,485,863       ~10 ms
p(10⁹)    =         22,801,763,489      ~100 ms     (the billionth prime)
p(10¹²)   =     29,996,224,275,833       ~14 s      (the trillionth prime)
p(10¹³)   =    323,780,508,946,331       ~85 s
p(10¹⁴)   =  3,475,385,758,524,527       ~7 min     (matches OEIS A006988)
p(2×10¹⁴) =  7,093,600,525,704,677      ~13 min     (the float64 frontier)
```

*(Node 22, single thread, pure JavaScript — run `npm run bench` yourself.
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

### 2. Count — exact π at the guess

**Primary engine — Lucy_Hedgehog's algorithm** (O(x^¾) time, O(√x) space):
a dynamic programme over the ≤ 2√x distinct values of ⌊x/k⌋, with the
recurrence `S(v) ← S(v) − [S(⌊v/p⌋) − π(p−1)]`. No factoring, no primality
tests. Implemented with blocked updates over typed arrays (`Uint32` for the
hot table) — measured ~10¹⁰ recurrence steps/second.

**Verification engine — Lagarias–Miller–Odlyzko (1985)**, the algorithm family
behind every prime-counting world record: π(x) = φ(x,a) + a − 1 − P₂(x,a),
with φ split into *ordinary leaves* (Möbius-weighted wheel counts) and
*special leaves* (answered by a segmented sieve with a Fenwick tree), and P₂
folded into one ascending segmented sweep. To our knowledge this is the
**first complete LMO implementation that runs inside a web page**.

The twist: in JavaScript, measured head-to-head, Lucy's lower constants beat
LMO(α=1) ~2× at every size below ~10¹⁷ — so Lucy computes your answer, and
LMO's job is *proof*: two unrelated algorithms, two unrelated bug surfaces,
agreeing digit-for-digit (`node cli.js --pi 1e13 --engine both`).

### 3. Walk — segmented sieve over the tiny gap

Exact π(guess) + the guess's tiny error ⇒ an odds-only segmented sieve of
Eratosthenes walks at most a few million integers to the exact n-th prime.
Milliseconds.

## How it is verified

`npm test` (~10 s) runs 10 independent layers, 172 checks (179 in --slow):

1. `primesUpTo` vs brute-force trial division;
2. table/sieve paths vs an independently generated prime list;
3. counting path vs sieve path across the cutover and 25 pseudo-random n;
4. exact π(10^k), k = 1…11, vs published values (k ≤ 14 in `--huge`);
5. p(10^k), k = 1…9, vs [OEIS A006988][oeis] (k ≤ 13 in `--slow`/`--huge`) —
   p(10⁹) and π(10¹²) were re-confirmed against independent web sources;
6. |R(x) − π(x)| ≪ √x (guards the Gram series / ζ implementation);
7. estimates always inside the rigorous Rosser/Dusart bracket;
8. input-validation edge cases;
9. **Lucy vs LMO cross-engine agreement** on structured + random x and
   several α values — the flagship check;
10. the self-contained `index.html` embeds the current engine byte-for-byte.

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

Type `1234567`, `1,234,567`, `1e9` or `10^12` and press **Find p(n)**.
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
| `cli.js` | nth-prime and π(x) command line |
| `test/test.js` | the verification suite |
| `bench.js` | timing table |

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
