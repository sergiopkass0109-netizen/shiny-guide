# The *n*-th Prime

Type a number **n** into the box, get back the **n-th prime number** — exactly,
for any **1 ≤ n ≤ 10¹²**, in milliseconds-to-seconds, with zero dependencies.

```
p(10⁶)  =          15,485,863      ~10 ms
p(10⁹)  =      22,801,763,489     ~120 ms      (the billionth prime)
p(10¹⁰) =     252,097,800,623     ~0.5 s
p(10¹¹) =   2,760,727,302,517     ~2.9 s
p(10¹²) =  29,996,224,275,833     ~16 s        (the trillionth prime)
```

*(Node 22, single thread, pure JavaScript. Run `npm run bench` yourself.)*

## Using it

**Web page** — `index.html` is one fully self-contained file: open it in any
browser (double-click works, no server needed), host it on GitHub Pages, or
serve it locally:

```sh
python3 -m http.server 8000      # or: npm run serve
# then visit http://localhost:8000
```

Enter `n` in the numerical box (plain `1234567`, `1,234,567`, `1e9` and
`10^12` all work), press **Find p(n)**, and the answer appears in the answer
box along with the diagnostics: the Riemann-R guess, the exact prime count at
the guess, how few primes the guess missed by, and the timing.

**Command line**

```sh
node cli.js 1e9          # p(1,000,000,000) = 22,801,763,489
node cli.js 10^12 --json
```

**Library**

```js
const { nthPrime, primeCount, riemannR } = require("./nthprime.js");
nthPrime(1e9).value   // 22801763489
primeCount(1e12)      // 37607912018  (exact π(10^12))
```

**Tests / benchmarks**

```sh
npm test               # 166-check verification suite (~4 s)
npm run test:slow      # adds p(10^10), p(10^11), π(10^12)  (~25 s)
npm run bench:big      # full timing table up to n = 10^12
```

## Why this is fast — the mathematics

Most "nth prime" code sieves all primes from 2 upward and counts. That is
O(p(n) log log p(n)) work and O(p(n)) memory — hopeless around n ≈ 10⁹ and
beyond. The fastest known approach (used by [primecount], which holds the
π(10²⁹) world record) never enumerates the primes at all. This project is a
faithful, dependency-free scale model of that architecture, in three acts:

### 1. Estimate — Cipolla seed, Riemann R refinement

The prime number theorem gives first-order estimates of p(n);
[Cipolla (1902)][cipolla] made them an asymptotic expansion (L = ln n,
LL = ln ln n):

```
p(n) ≈ n·( L + LL − 1 + (LL−2)/L − (LL² − 6·LL + 11)/(2L²) + … )
```

That seed is refined by Newton-inverting **Riemann's R function**, the
best classical smooth approximation to the prime-counting function π(x).
R(x) is evaluated with the [Gram series][gram], which converges for all
x > 0 (ζ computed via Euler–Maclaurin):

```
R(x) = 1 + Σ_{k≥1}  (ln x)^k / ( k · k! · ζ(k+1) )
```

Solving R(x) = n with Newton steps `x ← x − (R(x) − n)·ln x` gives a guess
x₀ for p(n) that is *astonishingly* good — measured on this implementation:

| n | guess missed the target by |
|---|---|
| 10⁹ | 1,429 primes (≈ 6·10⁻⁵ %) |
| 10¹² | 35,884 primes (≈ 4·10⁻⁶ %) |

Correctness never depends on the guess: x₀ is clamped into the **rigorous
bracket** `n(ln n + ln ln n − 1) ≤ p(n) < n(ln n + ln ln n)` (lower bound
[Dusart 1999][dusart], upper bound Rosser's theorem, n ≥ 6), and steps 2–3
locate the answer exactly wherever it is.

### 2. Count — Lucy_Hedgehog's exact π(x) in O(x¾) time, O(√x) space

The magic step: count primes ≤ x₀ **exactly without finding them**.
The method, posted by Project Euler user [Lucy_Hedgehog][lucy] in 2013, is a
dynamic programme over the ≤ 2√x distinct values of ⌊x/k⌋. Let S(v) be the
count of integers in [2, v] that survive sieving by every prime ≤ p; running
one Eratosthenes round with the next prime p removes exactly the composites
whose *least* prime factor is p, giving the recurrence

```
S(v) ← S(v) − [ S(⌊v/p⌋) − π(p−1) ]        for every v ≥ p²
```

Maintaining S only on the values {⌊x/k⌋} (two flat Float64 arrays: one
indexed by v ≤ √x, one by k = x/v ≤ √x) and applying the recurrence for each
prime p ≤ √x yields **π(x) = S(x)** in O(x^(3/4)) operations and O(√x)
memory — for x ≈ 3·10¹³ that is two 45 MB arrays instead of a 30-terabyte
sieve. Every quantity stays an exact integer below 2⁵³, so IEEE-754 doubles
are exact (the README header of `nthprime.js` carries the proof that
`Math.floor(a/b)` is exact for integer doubles a < 2⁵³).

This implementation computes exact values of π in:

| x | π(x), exact | time |
|---|---|---|
| 10¹⁰ | 455,052,511 | 73 ms |
| 10¹¹ | 4,118,054,813 | 0.26 s |
| 10¹² | 37,607,912,018 | 1.4 s |

*(For even larger x one would graduate to Meissel–Lehmer/LMO/Gourdon,
O(x^(2/3)/log²x) — the same idea with a cleverer recursion; see
[primecount]. At our 10¹² input cap, O(x¾) is the sweet spot of speed ×
provable simplicity.)*

### 3. Walk — segmented sieve of Eratosthenes over the tiny gap

Knowing exact π(x₀) and needing rank n, the answer lies |n − π(x₀)| primes
away — a few thousand, thanks to step 1. An **odds-only segmented sieve**
(windows sized adaptively from the known deficit, base primes ≤ √(bound) from
a simple sieve) walks to the exact prime. This step costs milliseconds.

### Numerical safety

* Inputs capped at n = 10¹² ⇒ p(n) < 3.1·10¹³ ≪ 2⁵³: every integer in the
  program is exactly representable as a double, and the floor-division
  exactness proof covers every quotient taken.
* No 32-bit bitwise tricks on any value that can exceed 2³¹.
* `isqrt`/`idiv` carry correction steps so a 1-ulp sqrt/division error can
  never change a result.

## How it is verified

`npm test` runs 166 checks (`--slow` adds the heavyweights) in 8 independent
layers:

1. `primesUpTo` against brute-force trial division;
2. the table/sieve paths against an independently generated prime list
   (exhaustive n = 1…2000, strided through n = 500,000);
3. **two independent algorithms against each other** — counting path vs
   sieve path — across the cutover boundary and 25 pseudo-random n;
4. exact π(10^k) for k = 1…12 against published values;
5. p(10^k) for k = 1…11 against [OEIS A006988][oeis] — the anchors
   p(10⁹) = 22,801,763,489 and π(10¹²) = 37,607,912,018 were re-confirmed
   against independent web sources for this project;
6. |R(x) − π(x)| stays ≪ √x (catches any Gram-series/ζ regression);
7. the estimate always lands inside the rigorous Rosser/Dusart bracket;
8. input-validation edge cases.

## Files

| file | purpose |
|---|---|
| `index.html` | the page: numerical box in, answer box out — **generated**, self-contained (computes in a Web Worker built from the inlined engine, with main-thread fallback) |
| `index.template.html` | source template for the page |
| `style.css` | styling (inlined into index.html) |
| `nthprime.js` | the engine (inlined into index.html; also used directly by Node) |
| `build.js` | `npm run build` → regenerates index.html from the three files above |
| `cli.js` | `node cli.js 1e9` |
| `test/test.js` | the verification suite (includes an index.html-in-sync check) |
| `bench.js` | timing table |

## Research sources

* [primecount — fast prime counting library][primecount] (architecture:
  R⁻¹ estimate + count + sieve; algorithms by Meissel, Lehmer,
  Lagarias–Miller–Odlyzko, Deléglise–Rivat, Gourdon)
* [Lucy_Hedgehog's O(x¾) prime counting][lucy] (Project Euler #10 thread;
  exposition by [Jacob Elafandi][lucy] and
  [griff's math blog](https://gbroxey.github.io/blog/2023/04/09/lucy-fenwick.html))
* [Gram series / Riemann prime counting function — MathWorld][gram]
* [Cipolla, "La determinazione assintotica dell'nᵐᵒ numero primo" (1902)][cipolla];
  modern treatment in [Axler, "New estimates for the nth prime number" (2019)](https://cs.uwaterloo.ca/journals/JIS/VOL22/Axler/axler17.pdf)
* [Dusart, "The k-th prime is greater than k(ln k + ln ln k − 1)" (1999)][dusart]
* [OEIS A006988 — the 10^n-th primes][oeis] and
  [How many primes are there? (Prime Pages)](https://t5k.org/howmany.html)

[primecount]: https://github.com/kimwalisch/primecount
[lucy]: https://math.berkeley.edu/~elafandi/euler/p10/
[gram]: https://mathworld.wolfram.com/GramSeries.html
[cipolla]: https://arxiv.org/pdf/1203.5413
[dusart]: https://t5k.org/howmany.html
[oeis]: https://oeis.org/A006988
