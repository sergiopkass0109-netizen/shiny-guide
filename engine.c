/* engine.c — Lucy_Hedgehog exact prime counting in freestanding C,
 * compiled to WebAssembly (see wasmbuild.js).  Same recurrence as the
 * verified primeCount() in nthprime.js; the test suite requires
 * digit-for-digit agreement between this, the JS Lucy engine and the JS LMO
 * engine on every value it checks.
 *
 * No libc.  Memory is the module's exported linear memory; the caller grows
 * it and passes byte offsets:
 *   small:   u32[r+1]     at smallOff    S(v) for v <= r,      r = isqrt(x)
 *   large:   u64[r+1]     at largeOff    S(floor(x/i)) for i <= r
 *   scratch: u64[SEG+2]   at scratchOff  per-segment difference array
 *
 * How the hot loop is organised (all of it exact; see the comments):
 *  - a pass for prime p updates large[i] -= S(x/(ip)) - pi(p-1); for
 *    i > r/p the value S(x/(ip)) is small[floor(xp/i)], xp = floor(x/p)
 *  - beyond sqrt(xp) consecutive i share a quotient q, and the run of q is
 *    (floor(xp/(q+1)), floor(xp/q)]: iterating q downward costs one
 *    independent division per run (no dependency chain)
 *  - primes are processed in blocks of up to KMAX consecutive primes
 *    p1 < ... < pk <= 2 p1.  For i > I0 = max(r/p1, x/p1^3) every read of
 *    every prime in the block is small[q] with q < p1^2 — a region no
 *    small-pass of the block touches — so the k passes are applied to
 *    large[(I0, imax]] in ONE cache-blocked sweep (segments of SEG entries),
 *    each prime's runs going to a difference array that one prefix-sum pass
 *    per segment materialises.  The prefix i <= I0 is then done prime by
 *    prime in the classical order (reads of large[m] with I0 < m <= r add
 *    back the sweep's contributions of the later primes).
 *
 * Build: clang --target=wasm32 -O3 -nostdlib -fno-builtin -mnontrapping-fptoint
 *        [-msimd128]                                        (see wasmbuild.js)
 */
#include <stdint.h>

typedef uint64_t u64;
typedef uint32_t u32;

#ifndef KMAX
#define KMAX 32      /* primes per block */
#endif
#ifndef SEG
#define SEG 4096     /* large[] entries per sweep segment (32 KiB) */
#endif

__attribute__((export_name("seg"))) u32 seg_export(void) { return SEG; }

/* small[0..len) -= sub — a run of equal updates (auto-vectorised in the SIMD build) */
static inline void sub_block(u32 *ptr, u64 len, u32 sub) {
  for (u64 k = 0; k < len; k++) ptr[k] -= sub;
}

static inline u64 isqrt64(u64 n) {
  if (n == 0) return 0;
  u64 r = (u64)__builtin_sqrt((double)n);
  while (r > 0 && r * r > n) r--;
  while ((r + 1) * (r + 1) <= n) r++;
  return r;
}

/* floor(a/b) via double division — exact for a < 2^53 (the IEEE proof in
 * nthprime.js applies verbatim; every caller stays below 9e15). */
static inline u64 fdiv(u64 a, u64 b) {
  return (u64)((double)a / (double)b);
}

__attribute__((export_name("isqrt")))
u64 isqrt_export(u64 n) { return isqrt64(n); }

/* progress callback into the host (fraction = p / r), every ~2048 primes */
__attribute__((import_name("progress"))) void host_progress(u64 p, u64 r);

/* large[i] -= small[floor(xp/i)] - sp1 for i in [a, b], where g = isqrt(xp).
 * i < g: one independent division per element.  i >= g: iterate the
 * quotient q downward; the run of q is (floor(xp/(q+1)), floor(xp/q)].
 * Exact: i <= floor(xp/q) => xp/i >= q, and i > floor(xp/(q+1)) =>
 * i > xp/(q+1) => xp/i < q+1, so floor(xp/i) == q on the whole run. */
static inline void tail_range(const u32 *small, u64 *large, u64 a, u64 b, u64 xp, u64 sp1, u64 g) {
  u64 i = a;
  u64 h = g; if (h < a) h = a; if (h > b + 1) h = b + 1;
  for (; i < h; i++) large[i] -= (u64)small[fdiv(xp, i)] - sp1;
  if (i > b) return;
  u64 q = fdiv(xp, i), qmin = fdiv(xp, b);
  if (qmin < 1) qmin = 1;
  u64 ePrev = i - 1;
  for (; q >= qmin; q--) {
    u64 e = fdiv(xp, q);
    if (e > b) e = b;
    u64 c = (u64)small[q] - sp1;
    for (u64 j = ePrev + 1; j <= e; j++) large[j] -= c;
    ePrev = e;
  }
}

/* Sweep version of tail_range: the runs go to the difference array D
 * (indexed i - base): D[start] += c, D[end+1] -= c.  A prefix sum over D
 * then equals, at every i, the total of all runs covering i — for every
 * prime of the block at once.  Head elements (i < g) are written directly. */
static inline void tail_range_d(const u32 *small, u64 *large, u64 a, u64 b, u64 xp, u64 sp1, u64 g, u64 *D, u64 base) {
  u64 i = a;
  u64 h = g; if (h < a) h = a; if (h > b + 1) h = b + 1;
  for (; i < h; i++) large[i] -= (u64)small[fdiv(xp, i)] - sp1;
  if (i > b) return;
  u64 q = fdiv(xp, i), qmin = fdiv(xp, b);
  if (qmin < 1) qmin = 1;
  u64 ePrev = i - 1;
  u64 *Db = D - base;
  for (; q >= qmin; q--) {
    u64 e = fdiv(xp, q);
    if (e > b) e = b;
    u64 c = (u64)small[q] - sp1;
    Db[ePrev + 1] += c;   /* u64 arithmetic wraps; the final sums are exact mod 2^64 */
    Db[e + 1] -= c;
    ePrev = e;
  }
}

/* small[v] -= small[floor(v/p)] - sp1 for v in [p^2, r], descending so that
 * small[q] (q = v/p < v) still holds the previous pass's value when read;
 * values of v sharing one quotient form a block. */
static inline void small_pass(u32 *small, u64 r, u64 p, u64 sp1) {
  u64 p2 = p * p;
  if (p2 > r) return;
  for (u64 v = r; v >= p2;) {
    u64 q = fdiv(v, p);
    u32 sub = (u32)((u64)small[q] - sp1);
    u64 w = q * p;
    if (w < p2) w = p2;
    sub_block(small + w, v - w + 1, sub);
    v = w - 1;
  }
}

__attribute__((export_name("pi_lucy")))
u64 pi_lucy(u64 x, u32 smallOff, u32 largeOff, u32 scratchOff) {
  if (x < 2) return 0;
  u64 r = isqrt64(x);
  u32 *small = (u32 *)(uintptr_t)smallOff;
  u64 *large = (u64 *)(uintptr_t)largeOff;
  u64 *D = (u64 *)(uintptr_t)scratchOff;

  for (u64 v = 1; v <= r; v++) small[v] = (u32)(v - 1);
  for (u64 i = 1; i <= r; i++) large[i] = fdiv(x, i) - 1;

  /* blocks start once the prefix I0 = max(r/p, x/p^3) is at most r/2,
   * i.e. p^3 >= 2x/r (a few hundred at most: the loop is negligible) */
  u64 pblock = 3, pb3 = 2 * (x / r);
  while (pblock * pblock * pblock < pb3) pblock++;

  u64 P[KMAX], XP[KMAX], SP1[KMAX], IMAX[KMAX], G[KMAX];
  u64 nextProg = 0;
  for (u64 p = 2; p <= r;) {
    if (p >= nextProg) { host_progress(p, r); nextProg = p + 2048; }
    if (small[p] == small[p - 1]) { p++; continue; }   /* p composite */

    if (p < pblock) {
      /* classical single-prime pass */
      u64 sp1 = small[p - 1], p2 = p * p, xp = fdiv(x, p);
      u64 imax = fdiv(x, p2); if (imax > r) imax = r;
      u64 isw = r / p; if (isw > imax) isw = imax;
      u64 i = 1;
      for (; i <= isw; i++) large[i] -= large[i * p] - sp1;
      tail_range(small, large, i, imax, xp, sp1, isqrt64(xp));
      small_pass(small, r, p, sp1);
      p++;
      continue;
    }

    /* ---- block of consecutive primes p = P[0] < ... < P[k-1] <= 2p ----
     * every P[j] - 1 and P[j] is below p^2, so small[] is final there and
     * both the primality test and SP1 = pi(P[j]-1) are already exact */
    u64 k = 0;
    for (u64 q = p; q <= r && q <= 2 * p && k < KMAX; q++)
      if (small[q] != small[q - 1]) P[k++] = q;
    u64 p1 = P[0];
    u64 I0 = r / p1;
    if ((double)p1 * (double)p1 * (double)p1 <= (double)x) {
      u64 t = fdiv(x, p1 * p1 * p1);
      if (t > I0) I0 = t;
    }
    for (u64 j = 0; j < k; j++) {
      XP[j] = fdiv(x, P[j]);
      IMAX[j] = fdiv(x, P[j] * P[j]); if (IMAX[j] > r) IMAX[j] = r;
      SP1[j] = small[P[j] - 1];
      G[j] = isqrt64(XP[j]);
    }

    /* 1. one cache-blocked sweep of large[(I0, IMAX[0]]] for all k primes.
     *    i > r/p1 => i*P[j] > r, so every read is small[floor(XP[j]/i)];
     *    i > x/p1^3 => that index is < p1^2 <= P[j]^2, a region untouched by
     *    every small-pass of this block — so all k passes commute here. */
    for (u64 L = I0 + 1; L <= IMAX[0];) {
      u64 R = L + SEG - 1; if (R > IMAX[0]) R = IMAX[0];
      u64 n = R - L + 1;
      for (u64 t = 0; t <= n; t++) D[t] = 0;
      for (u64 j = 0; j < k && IMAX[j] >= L; j++) {
        u64 Rj = R < IMAX[j] ? R : IMAX[j];
        tail_range_d(small, large, L, Rj, XP[j], SP1[j], G[j], D, L);
      }
      u64 acc = 0;
      for (u64 t = 0; t < n; t++) { acc += D[t]; large[L + t] -= acc; }
      L = R + 1;
    }

    /* 2. prefix i <= I0, prime by prime in the classical order, each
     *    followed by its small-pass (so small[] is in the right state for
     *    the next prime's reads at indices >= p1^2) */
    for (u64 j = 0; j < k; j++) {
      u64 pj = P[j], sp1 = SP1[j], xp = XP[j];
      u64 imaxj = IMAX[j] < I0 ? IMAX[j] : I0;
      u64 isw = r / pj; if (isw > imaxj) isw = imaxj;
      u64 i = 1;
      u64 i1 = I0 / pj; if (i1 > isw) i1 = isw;
      for (; i <= i1; i++) large[i] -= large[i * pj] - sp1;     /* i*pj <= I0: untouched by the sweep */
      for (; i <= isw; i++) {
        /* I0 < m <= r: the sweep already applied every prime l >= j whose
         * imax reaches m; adding those back gives S_{p_j - 1}(x/m) */
        u64 m = i * pj;
        u64 val = large[m];
        for (u64 l = j; l < k && IMAX[l] >= m; l++) val += (u64)small[fdiv(XP[l], m)] - SP1[l];
        large[i] -= val - sp1;
      }
      tail_range(small, large, i, imaxj, xp, sp1, G[j]);
      small_pass(small, r, pj, sp1);
    }
    p = P[k - 1] + 1;
  }
  return large[1];
}
