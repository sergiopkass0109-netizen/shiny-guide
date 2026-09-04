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
 *   wheel:   u16[30030]   at tblOff      cumulative coprime-to-30030 counts
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
 *    per segment materialises.  The prefix i <= I0 is done prime by prime in
 *    the classical order, half of the block before the sweep and half after,
 *    so that a read of large[m] with I0 < m <= r needs at most k/2
 *    corrections for the sweep's contributions of the other primes
 *
 * Build: clang --target=wasm32 -O3 -nostdlib -fno-builtin -mnontrapping-fptoint
 *        [-msimd128]                                        (see wasmbuild.js)
 */
#include <stdint.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t u8;
typedef int8_t i8;

/* Wheel start: the tables begin at the state after the passes for 2, 3, 5,
 * 7, 11 and 13 instead of running those six full passes.  Lucy's invariant
 * is S_p(v) = #{m in [2, v] : m prime or lpf(m) > p}; after 13 that is the
 * count of m <= v coprime to 30030, minus one for m = 1, plus the wheel
 * primes <= v.  T[m] = #{1 <= t <= m : gcd(t, 30030) = 1} for m < 30030. */
#define WHEEL 30030u
#define WHEEL_PHI 5760u
static void wheel_table(u16 *T) {
  u32 c = 0;
  T[0] = 0;
  for (u32 m = 1; m < WHEEL; m++) {
    if (m % 2 && m % 3 && m % 5 && m % 7 && m % 11 && m % 13) c++;
    T[m] = (u16)c;
  }
}
static inline u64 s13(u64 v, const u16 *T) {
  u64 pw = v >= 13 ? 6 : v >= 11 ? 5 : v >= 7 ? 4 : v >= 5 ? 3 : v >= 3 ? 2 : v >= 2 ? 1 : 0;
  return (v / WHEEL) * WHEEL_PHI + T[v % WHEEL] - 1 + pw;
}

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

/* prefix pass of block prime j over i in [1, imaxj], imaxj <= min(I0, IMAX[j]).
 * Reads of large[i*p]: direct while i*p <= I0; for I0 < i*p <= r the value
 * S_{p_j - 1}(x/m) is recovered from the sweep region — before the sweep by
 * subtracting the block's earlier primes, after it by adding back the later
 * ones (each T_l(m) = small[floor(XP_l/m)] - SP1_l, an index < p1^2, frozen). */
static inline void prefix_pass(const u32 *small, u64 *large, u64 r, u64 I0, u64 k, u64 j, int afterSweep,
                               const u64 *P, const u64 *XP, const u64 *SP1, const u64 *IMAX, const u64 *G, u64 imaxj) {
  u64 pj = P[j], sp1 = SP1[j], xp = XP[j];
  u64 isw = r / pj; if (isw > imaxj) isw = imaxj;
  u64 i = 1;
  u64 i1 = I0 / pj; if (i1 > isw) i1 = isw;
  for (; i <= i1; i++) large[i] -= large[i * pj] - sp1;     /* i*pj <= I0: untouched by the sweep */
  for (; i <= isw; i++) {
    u64 m = i * pj;
    u64 val = large[m];
    if (afterSweep) { for (u64 l = j; l < k && IMAX[l] >= m; l++) val += (u64)small[fdiv(XP[l], m)] - SP1[l]; }
    else            { for (u64 l = 0; l < j && IMAX[l] >= m; l++) val -= (u64)small[fdiv(XP[l], m)] - SP1[l]; }
    large[i] -= val - sp1;
  }
  tail_range(small, large, i, imaxj, xp, sp1, G[j]);
}

__attribute__((export_name("pi_lucy")))
u64 pi_lucy(u64 x, u32 smallOff, u32 largeOff, u32 scratchOff, u32 tblOff) {
  if (x < 2) return 0;
  u64 r = isqrt64(x);
  u32 *small = (u32 *)(uintptr_t)smallOff;
  u64 *large = (u64 *)(uintptr_t)largeOff;
  u64 *D = (u64 *)(uintptr_t)scratchOff;
  u16 *T = (u16 *)(uintptr_t)tblOff;

  /* tables at the state after the primes <= 13; for v < 17^2 that is pi(v),
   * so the primality tests below are already exact from p = 17 on */
  wheel_table(T);
  for (u64 v = 1; v <= r; v++) small[v] = (u32)s13(v, T);
  for (u64 i = 1; i <= r; i++) large[i] = s13(fdiv(x, i), T);

  /* blocks start once the prefix I0 = max(r/p, x/p^3) is at most r/2,
   * i.e. p^3 >= 2x/r (a few hundred at most: the loop is negligible) */
  u64 pblock = 3, pb3 = 2 * (x / r);
  while (pblock * pblock * pblock < pb3) pblock++;

  u64 P[KMAX], XP[KMAX], SP1[KMAX], IMAX[KMAX], G[KMAX];
  u64 nextProg = 0;
  for (u64 p = 17; p <= r;) {
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

    u64 half = k / 2;

    /* 1. prefixes of the first half, in the classical order, each followed
     *    by its small-pass (so small[] is in the right state for the next
     *    prime's reads at indices >= p1^2); the sweep has not run yet */
    for (u64 j = 0; j < half; j++) {
      prefix_pass(small, large, r, I0, k, j, 0, P, XP, SP1, IMAX, G, IMAX[j] < I0 ? IMAX[j] : I0);
      small_pass(small, r, P[j], SP1[j]);
    }

    /* 2. one cache-blocked sweep of large[(I0, IMAX[0]]] for all k primes.
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

    /* 3. prefixes of the second half, now correcting for the sweep */
    for (u64 j = half; j < k; j++) {
      prefix_pass(small, large, r, I0, k, j, 1, P, XP, SP1, IMAX, G, IMAX[j] < I0 ? IMAX[j] : I0);
      small_pass(small, r, P[j], SP1[j]);
    }
    p = P[k - 1] + 1;
  }
  return large[1];
}

/* ======================================================================
 * Deléglise–Rivat prime counting (the Lagarias–Miller–Odlyzko family):
 * O(x^(2/3) / log^2 x) time, O(sqrt(x)) memory.  Exact.
 *
 *   pi(x) = S0 + S_special + a - 1 - P2(x, a),   y = alpha * cbrt(x), a = pi(y)
 *   S0        = sum_{n <= y, n = 1 or (mu(n) != 0, lpf(n) > 13)} mu(n) phi6(x/n)
 *   S_special = sum over leaves (n, b): n <= y < n p_b, mu(n) != 0, 13 < p_b < lpf(n),
 *               b <= a, of  -mu(n) * phi(floor(x / (n p_b)), b - 1)
 *   P2        = sum_{y < p <= sqrt(x)} (pi(x/p) - pi(p) + 1)
 *
 * Special leaves split by v = floor(x / (n p_b)):
 *   v <  p_b    trivial   phi = 1                    (closed form; p_b > sqrt(y) => n prime)
 *   v <  p_b^2  easy      phi = pi(v) - b + 2         (pi from the fully sieved segment)
 *   otherwise   hard      phi = count of numbers <= v unsieved by p_1..p_{b-1}
 * For p_b <= sqrt(y) every leaf is answered from the sieve state (exact for any v).
 *
 * One segmented bit sieve over [1, z], z = x/y: within a segment the primes are
 * crossed off in order (the primes <= 149 word-wise through residue masks) and
 * each prime's hard leaves are answered just before it is crossed off, from
 * per-block popcount counters walked monotonically; once the segment is fully
 * sieved by the primes <= sqrt(z), its unsieved numbers are 1 and the primes
 * above sqrt(z), which answers the easy leaves and P2 in O(1) each.
 *
 * Memory is carved from `base`: dr_bytes(x, y) tells the host how much.  The
 * host picks y (see drAlpha in nthprime.js); 17 <= y < sqrt(x) and a >= 8 are
 * required, otherwise the call returns ~0 and the host uses the Lucy engine.
 * Verified digit-for-digit against the Lucy engine (see test/test.js).
 * ====================================================================== */
#define DRS (1u << 20)              /* numbers per segment */
#define DRW (DRS / 64)              /* words per segment */
#define DRNB (DRW / 64)             /* big blocks (4096 numbers) per segment */
#define P_C 13                      /* wheel primes 2,3,5,7,11,13: c = 6 */
#define PMASK 149                   /* primes 17..PMASK are crossed off word-wise with residue masks */

static inline u64 umin(u64 a, u64 b) { return a < b ? a : b; }
static inline u64 umax(u64 a, u64 b) { return a > b ? a : b; }
static inline u64 align8(u64 v) { return (v + 7) & ~(u64)7; }

static u64 G[8];                                        /* components of the last call (dr_part), for the test suite */
__attribute__((export_name("dr_part"))) u64 dr_part(u32 k) { return G[k]; }

/* upper bound on pi(n) for the memory budget: 1.26 n / ln n + 64 (n >= 17) */
static u64 pi_bound(u64 n) {
  if (n < 17) return 80;
  u64 lg2 = 63 - (u64)__builtin_clzll(n);              /* floor(log2 n) <= log2 n, so this only over-allocates */
  double ln = 0.6931471805599453 * (double)lg2;
  return (u64)(1.26 * (double)n / ln) + 64;
}

/* bytes needed by pi_dr for (x, y); the host grows memory to this */
__attribute__((export_name("dr_bytes")))
u64 dr_bytes(u64 x, u64 y) {
  u64 r = isqrt64(x);
  u64 nb = pi_bound(r);
  u64 bytes = 0;
  bytes += align8(4 * (nb + 2));           /* prime list */
  bytes += align8(8 * (r / 64 + 3));       /* prime bit sieve up to r */
  bytes += align8(y + 2);                  /* mu */
  bytes += align8(4 * (y + 2));            /* lpf */
  bytes += align8(4 * (y + 2));            /* pi */
  bytes += align8(2 * WHEEL);              /* wheel table */
  bytes += align8(8 * 64);                 /* wheel masks */
  bytes += align8(8 * 4096);               /* residue masks for the primes 17..PMASK (sum of p = 2235) */
  bytes += align8(8 * 64);                 /* mask offsets per prime index */
  bytes += align8(8 * DRW);                /* sieve words */
  bytes += align8(2 * DRNB);               /* cbig */
  bytes += align8(4 * (DRW + 1));          /* cum */
  u64 a = pi_bound(y);
  bytes += align8(8 * (a + 2)) * 4;        /* cnt, nxt, nS, qIdx */
  return bytes + 4096;
}

/* Running count of unsieved numbers at offsets 0..pos of the segment, for a
 * monotone sequence of pos (the leaves of one prime arrive in increasing
 * order): the sum of whole big blocks (cbig, maintained during cross-off)
 * plus popcounts of the words of the current block scanned once. */
typedef struct { u64 B, sb, w, sw; } RunCount;
static inline void rc_init(RunCount *rc) { rc->B = 0; rc->sb = 0; rc->w = 0; rc->sw = 0; }
static inline u64 rc_count(RunCount *rc, const u64 *sv, const u16 *cbig, u64 pos) {
  u64 nB = pos >> 12, nw = pos >> 6;
  if (nB != rc->B) { while (rc->B < nB) rc->sb += cbig[rc->B++]; rc->w = nB << 6; rc->sw = 0; }
  while (rc->w < nw) rc->sw += (u64)__builtin_popcountll(sv[rc->w++]);
  u64 sh = pos & 63;
  u64 mask = sh == 63 ? ~(u64)0 : (((u64)1 << (sh + 1)) - 1);
  return rc->sb + rc->sw + (u64)__builtin_popcountll(sv[nw] & mask);
}

__attribute__((export_name("pi_dr")))
u64 pi_dr(u64 x, u64 y, u32 base) {
  if (x < 2) return 0;
  u64 r = isqrt64(x);
  if (y < 17 || y >= r) return ~(u64)0;                /* caller keeps 17 <= y < sqrt(x) */
  u64 z = fdiv(x, y);
  u64 off = base;
  u64 nb = pi_bound(r);
  u32 *PR = (u32 *)(uintptr_t)off; off += align8(4 * (nb + 2));
  u64 *PS = (u64 *)(uintptr_t)off; off += align8(8 * (r / 64 + 3));
  i8  *mu = (i8 *)(uintptr_t)off;  off += align8(y + 2);
  u32 *lpf = (u32 *)(uintptr_t)off; off += align8(4 * (y + 2));
  u32 *pi = (u32 *)(uintptr_t)off;  off += align8(4 * (y + 2));
  u16 *T = (u16 *)(uintptr_t)off;   off += align8(2 * WHEEL);
  u64 *MSK = (u64 *)(uintptr_t)off; off += align8(8 * 64);
  u64 *PMK = (u64 *)(uintptr_t)off; off += align8(8 * 4096);
  u64 *PMO = (u64 *)(uintptr_t)off; off += align8(8 * 64);
  u64 *sv = (u64 *)(uintptr_t)off;  off += align8(8 * DRW);
  u16 *cbig = (u16 *)(uintptr_t)off; off += align8(2 * DRNB);
  u32 *cum = (u32 *)(uintptr_t)off; off += align8(4 * (DRW + 1));
  u64 abound = pi_bound(y);
  u64 *cnt = (u64 *)(uintptr_t)off; off += align8(8 * (abound + 2));
  u64 *nxt = (u64 *)(uintptr_t)off; off += align8(8 * (abound + 2));
  u64 *nS = (u64 *)(uintptr_t)off;  off += align8(8 * (abound + 2));
  u64 *qIdx = (u64 *)(uintptr_t)off; off += align8(8 * (abound + 2));

  /* ---- primes up to r (bit sieve, 1 bit per number) ---- */
  u64 nw = r / 64 + 2;
  for (u64 w = 0; w < nw; w++) PS[w] = ~(u64)0;
  PS[0] &= ~(u64)3;                                     /* 0 and 1 */
  for (u64 i = 2; i * i <= r; i++)
    if ((PS[i >> 6] >> (i & 63)) & 1)
      for (u64 m = i * i; m <= r; m += i) PS[m >> 6] &= ~((u64)1 << (m & 63));
  u64 nP = 0;
  for (u64 i = 2; i <= r; i++) if ((PS[i >> 6] >> (i & 63)) & 1) { if (nP >= nb) return ~(u64)0; PR[nP++] = (u32)i; }

  /* ---- mu, lpf, pi up to y ---- */
  for (u64 n = 0; n <= y; n++) { mu[n] = 1; lpf[n] = 0; }
  for (u64 k = 0; k < nP && PR[k] <= y; k++) {
    u64 p = PR[k];
    for (u64 m = p; m <= y; m += p) { mu[m] = (i8)-mu[m]; if (lpf[m] == 0) lpf[m] = (u32)p; }
    for (u64 m = p * p; m <= y; m += p * p) mu[m] = 0;
  }
  pi[0] = 0;
  for (u64 n = 1; n <= y; n++) pi[n] = pi[n - 1] + (n >= 2 && lpf[n] == n);
  u64 a = pi[y];                                        /* 1-based count: p_a = largest prime <= y */
  if (a < 8) return ~(u64)0;                            /* need p_a > 13 */

  wheel_table(T);                                       /* T[m] = #{1 <= t <= m : gcd(t, 30030) = 1} */
  #define PHI6(v) (((v) / WHEEL) * (u64)WHEEL_PHI + T[(v) % WHEEL])

  /* ---- S0: ordinary leaves ---- */
  u64 S0 = 0;                                           /* modular u64: exact at the end */
  for (u64 n = 1; n <= y; n++) {
    if (n == 1) { S0 += PHI6(x); continue; }
    if (mu[n] == 0 || lpf[n] <= P_C) continue;
    u64 v = PHI6(fdiv(x, n));
    if (mu[n] > 0) S0 += v; else S0 -= v;
  }

  /* ---- index bounds ---- */
  u64 sqy = isqrt64(y);
  u64 bSmall = pi[sqy];                                 /* leaves of b <= bSmall: any n, via the sieve */
  if (bSmall < 6) bSmall = 6;                           /* special leaves start at b = 7 (p_7 = 17 > 13) */
  u64 x14 = isqrt64(r);                                 /* floor(x^(1/4)) */
  u64 bHard = umax(bSmall, pi[umin(x14, y)]);           /* hard leaves exist only for b <= bHard */
  u64 sqz = isqrt64(z);
  u64 aS = pi[umin(sqz, y)];                            /* primes crossed off in the segment sieve */
  if (aS < bHard) aS = bHard;

  /* ---- trivial leaves, closed form, for p_b > sqrt(y) ---- */
  u64 Striv = 0;
  for (u64 b = bSmall + 1; b <= a; b++) {
    u64 p = PR[b - 1];
    u64 lo = umax(umax(p, fdiv(y, p)), fdiv(fdiv(x, p), p));   /* q > lo */
    if (lo < y) Striv += pi[y] - pi[lo];
  }

  /* ---- masks for the wheel primes 3..13 (per residue of the word's first number) ---- */
  { static const u64 wp[5] = {3, 5, 7, 11, 13}; u64 o = 0;
    for (u64 k = 0; k < 5; k++) { u64 p = wp[k];
      for (u64 res = 0; res < p; res++) { u64 m = 0; for (u64 j = 0; j < 64; j++) if ((res + j) % p == 0) m |= (u64)1 << j; MSK[o + res] = m; }
      o += p; } }
  #define EVEN 0xAAAAAAAAAAAAAAAAull                    /* lo is odd: even numbers sit at odd offsets */
  /* residue masks for the primes 17..PMASK: PMK[PMO[b] + res] has the bits of the
   * multiples of p_b in a word whose first number is res (mod p_b) */
  u64 bMask = 6;
  { u64 o = 0;
    for (u64 b = 7; b <= aS && PR[b - 1] <= PMASK; b++) { u64 p = PR[b - 1]; PMO[b] = o;
      for (u64 res = 0; res < p; res++) { u64 m = 0; for (u64 j = 0; j < 64; j++) if ((res + j) % p == 0) m |= (u64)1 << j; PMK[o + res] = m; }
      o += p; bMask = b; } }

  /* ---- cursors ---- */
  for (u64 b = 0; b <= aS + 1; b++) { cnt[b] = 0; nxt[b] = 0; nS[b] = 0; qIdx[b] = 0; }
  for (u64 b = 7; b <= aS; b++) nxt[b] = PR[b - 1];      /* first multiple to cross off: p itself */
  for (u64 b = 7; b <= bSmall; b++) nS[b] = y;           /* small b: n descends from y */
  for (u64 b = bSmall + 1; b <= bHard; b++) {            /* large b: q descends from min(y, x/p^3) */
    u64 p = PR[b - 1];
    u64 qmax = umin(y, fdiv(fdiv(fdiv(x, p), p), p));
    qIdx[b] = qmax >= 2 ? pi[qmax] : 0;                  /* 1-based index of the largest prime <= qmax; 0 = none */
  }
  u64 p2Idx = nP;                                        /* P2 cursor: primes descending from the largest <= r (1-based index) */
  u64 Shard = 0, Seasy = 0, P2 = 0;
  u64 belowAll = 0;                                      /* unsieved numbers (fully sieved state) below the segment */
  u64 nseg = (z - 1) / DRS + 1;

  for (u64 seg = 0; seg < nseg; seg++) {
    u64 lo = 1 + seg * DRS;
    u64 hi = umin(lo + DRS - 1, z);
    u64 len = hi - lo + 1;
    if ((seg & 15) == 0) host_progress(seg, nseg);
    /* wheel state: multiples of 2..13 removed, 1 kept */
    { u64 r3 = lo % 3, r5 = lo % 5, r7 = lo % 7, r11 = lo % 11, r13 = lo % 13;
      for (u64 w = 0; w < DRW; w++) {
        u64 m = EVEN | MSK[r3] | MSK[3 + r5] | MSK[8 + r7] | MSK[15 + r11] | MSK[26 + r13];
        sv[w] = ~m;
        r3 = (r3 + 64) % 3; r5 = (r5 + 64) % 5; r7 = (r7 + 64) % 7; r11 = (r11 + 64) % 11; r13 = (r13 + 64) % 13;
      }
      if (len < DRS) {                                   /* clear offsets >= len */
        u64 fw = len >> 6, fb = len & 63;
        if (fb) sv[fw] &= ((u64)1 << fb) - 1, fw++;
        for (u64 w = fw; w < DRW; w++) sv[w] = 0;
      } }
    u64 curUnm = 0;
    for (u64 B = 0; B < DRNB; B++) { u64 s = 0; for (u64 t = 0; t < 64; t++) s += (u64)__builtin_popcountll(sv[(B << 6) + t]); cbig[B] = (u16)s; curUnm += s; }

    for (u64 b = 7; b <= aS; b++) {
      u64 p = PR[b - 1];
      if (b <= bHard) {
        /* hard leaves of b: state b-1 = cnt[b-1] (below lo) + count in segment */
        RunCount rc; rc_init(&rc);
        if (b <= bSmall) {
          u64 n = nS[b], nmin = fdiv(y, p);              /* leaves need n > y/p */
          while (n > nmin) {
            u64 v = fdiv(x, n * p);
            if (v > hi) break;
            if (mu[n] != 0 && lpf[n] > p) {
              u64 phi = cnt[b - 1] + rc_count(&rc, sv, cbig, v - lo);
              if (mu[n] > 0) Shard -= phi; else Shard += phi;
            }
            n--;
          }
          nS[b] = n;
        } else {
          u64 qi = qIdx[b], qmin = umax(p, fdiv(y, p));   /* leaves need q > qmin */
          while (qi >= 1 && PR[qi - 1] > qmin) {
            u64 q = PR[qi - 1];
            u64 v = fdiv(x, q * p);
            if (v > hi) break;
            Shard += cnt[b - 1] + rc_count(&rc, sv, cbig, v - lo);   /* -mu(q) = +1 */
            qi--;
          }
          qIdx[b] = qi;
        }
        cnt[b - 1] += curUnm;
        if (b <= bMask) {
          /* cross off p word-wise with residue masks, counted per word */
          const u64 *M = PMK + PMO[b];
          u64 res = lo % p, step = 64 % p, wEnd = (len + 63) >> 6;
          for (u64 w = 0; w < wEnd; w++) {                 /* branchless: every word is processed */
            u64 old = sv[w], hit = old & M[res];
            sv[w] = old & ~hit;
            u64 d = (u64)__builtin_popcountll(hit);
            cbig[w >> 6] -= (u16)d; curUnm -= d;
            res += step; res -= p & (0 - (u64)(res >= p));
          }
        } else {
          /* cross off p, counted */
          u64 m = nxt[b];
          for (; m <= hi; m += p) {
            u64 j = m - lo, w = j >> 6, sh = j & 63;
            u64 old = sv[w], d = (old >> sh) & 1;
            sv[w] = old & ~((u64)1 << sh);
            cbig[w >> 6] -= (u16)d; curUnm -= d;
          }
          nxt[b] = m;
        }
      } else if (b <= bMask) {
        const u64 *M = PMK + PMO[b];
        u64 res = lo % p, step = 64 % p, wEnd = (len + 63) >> 6;
        for (u64 w = 0; w < wEnd; w++) { sv[w] &= ~M[res]; res += step; res -= p & (0 - (u64)(res >= p)); }
      } else {
        /* cross off p, plain */
        u64 m = nxt[b];
        for (; m <= hi; m += p) { u64 j = m - lo; sv[j >> 6] &= ~((u64)1 << (j & 63)); }
        nxt[b] = m;
      }
    }

    /* fully sieved: cumulative counts; unsieved = 1 (first segment) and primes > p_aS */
    { u64 s = 0; for (u64 w = 0; w < DRW; w++) { cum[w] = (u32)s; s += (u64)__builtin_popcountll(sv[w]); } cum[DRW] = (u32)s; }
    #define UNS_LE(pos) (belowAll + cum[(pos) >> 6] + (u64)__builtin_popcountll(sv[(pos) >> 6] & (((pos) & 63) == 63 ? ~(u64)0 : (((u64)1 << (((pos) & 63) + 1)) - 1))))
    #define PI_AT(v) ((v) <= y ? (u64)pi[(v)] : aS + UNS_LE((v) - lo) - 1)

    /* P2: primes p in (y, r] with x/p in [lo, hi]; p descends as segments ascend */
    while (p2Idx > a) {
      u64 p = PR[p2Idx - 1];
      u64 v = fdiv(x, p);
      if (v > hi) break;
      P2 += PI_AT(v) - p2Idx + 1;                        /* pi(x/p) - pi(p) + 1, pi(p) = p2Idx */
      p2Idx--;
    }

    /* easy leaves: p_b > sqrt(y), n = q prime, p_b <= v < p_b^2, v in [lo, hi] */
    for (u64 b = bSmall + 1; b <= a; b++) {
      u64 p = PR[b - 1];
      u64 xp = fdiv(x, p), xp2 = fdiv(xp, p), xp3 = fdiv(xp2, p);
      u64 qlo = umax(umax(p, fdiv(y, p)), xp3) + 1;      /* q > p, q > y/p, q > x/p^3 */
      u64 qhi = umin(y, xp2);                            /* q <= x/p^2 (else trivial) */
      /* segment: lo <= floor(xp/q) <= hi  <=>  xp/(hi+1) < q <= xp/lo */
      u64 sl = fdiv(xp, hi + 1) + 1, sh2 = fdiv(xp, lo);
      if (sl > qlo) qlo = sl;
      if (sh2 < qhi) qhi = sh2;
      if (qlo > qhi) continue;
      u64 i0 = pi[qlo - 1] + 1, i1 = pi[qhi];             /* 1-based prime indices in [qlo, qhi] */
      for (u64 i = i0; i <= i1; i++) {
        u64 v = fdiv(xp, PR[i - 1]);
        Seasy += PI_AT(v) - b + 2;
      }
    }
    belowAll += cum[DRW];
  }
  host_progress(nseg, nseg);
  G[0] = S0; G[1] = Shard; G[2] = Seasy; G[3] = Striv; G[4] = P2; G[5] = a; G[6] = aS; G[7] = bHard;
  return S0 + Shard + Seasy + Striv + a - 1 - P2;
}
