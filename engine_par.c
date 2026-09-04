/* engine_par.c — range kernels for the multi-core engine (parallel.js).
 * The loops of engine.c restricted to an index range, so K threads can each
 * run their slice over ONE shared linear memory.  Built with atomics + shared
 * memory (see wasmbuild.js); the barriers and the schedule live in JS.
 *
 * The `large` table is stored as IEEE doubles (exact below 2^53, and every
 * value here is <= x < 9e15) so that the JavaScript coordinator can read and
 * write the same bytes through a Float64Array; `small` is u32 on both sides.
 * Mixing u64 and f64 views of one buffer was a bug we caught.
 *
 * Block parameters (doubles at bpOff, written by the coordinator):
 *   bp[0] = k (primes in the block)   bp[1] = I0 (prefix bound)   bp[2] = r
 *   bp[8 + 5j + 0..4] = P[j], XP[j] = floor(x/P[j]), SP1[j] = pi(P[j]-1),
 *                      IMAX[j] = min(r, floor(x/P[j]^2)), G[j] = isqrt(XP[j])
 * A classical single-prime pass is a block with k = 1 and I0 = r.
 *
 * Work distribution: the *_dyn kernels let every thread claim SEG-sized
 * chunks of a range through one atomic counter (ctrl[2]) until the range is
 * exhausted, so the division-heavy low end of a range never lands on a
 * single thread.  Each chunk is independent: the sweep has no hazards at
 * all, and a prefix tail only reads indices above the whole range.
 */
#include <stdint.h>

typedef uint64_t u64;
typedef uint32_t u32;

#ifndef SEG
#define SEG 4096     /* large[] entries per sweep segment (32 KiB); scratch is SEG+2 doubles per thread */
#endif

__attribute__((export_name("seg"))) u32 seg_export(void) { return SEG; }

static inline void sub_block(u32 *ptr, u64 len, u32 sub) {
  for (u64 k = 0; k < len; k++) ptr[k] -= sub;
}

/* exact floor(a/b) for a < 2^53 via double division (see engine.c) */
static inline u64 fdiv(u64 a, u64 b) {
  return (u64)((double)a / (double)b);
}

#define BP_P(j)    ((u64)bp[8 + 5 * (j)])
#define BP_XP(j)   ((u64)bp[9 + 5 * (j)])
#define BP_SP1(j)  ((u64)bp[10 + 5 * (j)])
#define BP_IMAX(j) ((u64)bp[11 + 5 * (j)])
#define BP_G(j)    ((u64)bp[12 + 5 * (j)])

/* large[i] -= small[floor(xp/i)] - s for i in [a, b]; g = isqrt(xp).
 * Head (i < g): one independent division per element; then runs of equal
 * quotient q, iterated downward — the run of q is (floor(xp/(q+1)), floor(xp/q)]. */
static inline void tail_range(const u32 *small, double *large, u64 a, u64 b, u64 xp, double s, u64 g) {
  u64 i = a;
  u64 h = g; if (h < a) h = a; if (h > b + 1) h = b + 1;
  for (; i < h; i++) large[i] -= (double)small[fdiv(xp, i)] - s;
  if (i > b) return;
  u64 q = fdiv(xp, i), qmin = fdiv(xp, b);
  if (qmin < 1) qmin = 1;
  u64 ePrev = i - 1;
  for (; q >= qmin; q--) {
    u64 e = fdiv(xp, q);
    if (e > b) e = b;
    double c = (double)small[q] - s;
    for (u64 j = ePrev + 1; j <= e; j++) large[j] -= c;
    ePrev = e;
  }
}

/* sweep version: runs go to the difference array D (indexed i - base) */
static inline void tail_range_d(const u32 *small, double *large, u64 a, u64 b, u64 xp, double s, u64 g, double *D, u64 base) {
  u64 i = a;
  u64 h = g; if (h < a) h = a; if (h > b + 1) h = b + 1;
  for (; i < h; i++) large[i] -= (double)small[fdiv(xp, i)] - s;
  if (i > b) return;
  u64 q = fdiv(xp, i), qmin = fdiv(xp, b);
  if (qmin < 1) qmin = 1;
  u64 ePrev = i - 1;
  double *Db = D - base;
  for (; q >= qmin; q--) {
    u64 e = fdiv(xp, q);
    if (e > b) e = b;
    double c = (double)small[q] - s;
    Db[ePrev + 1] += c;   /* partial sums stay integers far below 2^53: exact */
    Db[e + 1] -= c;
    ePrev = e;
  }
}

__attribute__((export_name("init_range")))
void init_range(u32 smallOff, u32 largeOff, u64 x, u64 a, u64 b) {
  u32 *small = (u32 *)(uintptr_t)smallOff;
  double *large = (double *)(uintptr_t)largeOff;
  for (u64 v = a; v <= b; v++) {
    small[v] = (u32)(v - 1);
    large[v] = (double)(fdiv(x, v) - 1);
  }
}

/* prefix pass of block prime j over i in [a, b] (b <= min(I0, IMAX[j])):
 * reads of large[i*p] are direct while i*p <= I0, corrected for the sweep's
 * contributions of primes l >= j while i*p <= r, and small[] beyond. */
__attribute__((export_name("prefix_range")))
void prefix_range(u32 smallOff, u32 largeOff, u32 bpOff, u32 j, u64 a, u64 b) {
  u32 *small = (u32 *)(uintptr_t)smallOff;
  double *large = (double *)(uintptr_t)largeOff;
  const double *bp = (const double *)(uintptr_t)bpOff;
  u64 k = (u64)bp[0], I0 = (u64)bp[1], r = (u64)bp[2];
  u64 p = BP_P(j), xp = BP_XP(j), g = BP_G(j);
  double s = (double)BP_SP1(j);
  u64 isw = r / p; if (isw > b) isw = b;
  u64 i1 = I0 / p; if (i1 > isw) i1 = isw;
  u64 i = a;
  for (; i <= i1; i++) large[i] -= large[i * p] - s;
  for (; i <= isw; i++) {
    u64 m = i * p;
    double val = large[m];
    for (u64 l = j; l < k && BP_IMAX(l) >= m; l++) val += (double)small[fdiv(BP_XP(l), m)] - (double)BP_SP1(l);
    large[i] -= val - s;
  }
  tail_range(small, large, i, b, xp, s, g);
}

/* cache-blocked sweep of large[[a, b]] (a > I0, b <= IMAX[0]) for every
 * prime of the block; D is this thread's scratch of SEG+2 doubles */
__attribute__((export_name("sweep_range")))
void sweep_range(u32 smallOff, u32 largeOff, u32 bpOff, u32 dOff, u64 a, u64 b) {
  u32 *small = (u32 *)(uintptr_t)smallOff;
  double *large = (double *)(uintptr_t)largeOff;
  const double *bp = (const double *)(uintptr_t)bpOff;
  double *D = (double *)(uintptr_t)dOff;
  u64 k = (u64)bp[0];
  for (u64 L = a; L <= b;) {
    u64 R = L + SEG - 1; if (R > b) R = b;
    u64 n = R - L + 1;
    for (u64 t = 0; t <= n; t++) D[t] = 0.0;
    for (u64 j = 0; j < k && BP_IMAX(j) >= L; j++) {
      u64 Rj = R < BP_IMAX(j) ? R : BP_IMAX(j);
      tail_range_d(small, large, L, Rj, BP_XP(j), (double)BP_SP1(j), BP_G(j), D, L);
    }
    double acc = 0.0;
    for (u64 t = 0; t < n; t++) { acc += D[t]; large[L + t] -= acc; }
    L = R + 1;
  }
}

/* small[v] -= small[floor(v/p)] - sp1 for v in [p^2, r] (see small_range) */
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

/* the whole step 2 of a block on one thread: prefix pass of every prime over
 * i <= min(I0, IMAX[j]) followed by its small-pass — used when the block's
 * prefixes are too small to be worth a barrier (one call instead of k) */
__attribute__((export_name("block_tail")))
void block_tail(u32 smallOff, u32 largeOff, u32 bpOff) {
  u32 *small = (u32 *)(uintptr_t)smallOff;
  const double *bp = (const double *)(uintptr_t)bpOff;
  u64 k = (u64)bp[0], I0 = (u64)bp[1], r = (u64)bp[2];
  for (u64 j = 0; j < k; j++) {
    u64 imaxj = BP_IMAX(j) < I0 ? BP_IMAX(j) : I0;
    if (imaxj >= 1) prefix_range(smallOff, largeOff, bpOff, (u32)j, 1, imaxj);
    small_pass(small, r, BP_P(j), BP_SP1(j));
  }
}

static inline u64 claim(u32 ctrlOff) {
  u32 *ctrl = (u32 *)(uintptr_t)ctrlOff;
  return (u64)__atomic_fetch_add(&ctrl[2], 1u, __ATOMIC_SEQ_CST);
}

/* sweep of [lo, hi] shared dynamically: chunk s is [lo + s*chunk, lo + (s+1)*chunk - 1] */
__attribute__((export_name("sweep_dyn")))
void sweep_dyn(u32 smallOff, u32 largeOff, u32 bpOff, u32 dOff, u32 ctrlOff, u64 lo, u64 hi, u64 chunk) {
  if (chunk < 1) chunk = 1;
  if (chunk > SEG) chunk = SEG;
  for (;;) {
    u64 L = lo + claim(ctrlOff) * chunk;
    if (L > hi) break;
    u64 R = L + chunk - 1; if (R > hi) R = hi;
    sweep_range(smallOff, largeOff, bpOff, dOff, L, R);
  }
}

/* hazard-free prefix tail of prime j over [lo, hi], shared dynamically */
__attribute__((export_name("prefix_dyn")))
void prefix_dyn(u32 smallOff, u32 largeOff, u32 bpOff, u32 ctrlOff, u32 j, u64 lo, u64 hi, u64 chunk) {
  if (chunk < 1) chunk = 1;
  for (;;) {
    u64 L = lo + claim(ctrlOff) * chunk;
    if (L > hi) break;
    u64 R = L + chunk - 1; if (R > hi) R = hi;
    prefix_range(smallOff, largeOff, bpOff, j, L, R);
  }
}

__attribute__((export_name("small_range")))
void small_range(u32 smallOff, u64 a, u64 b, u64 p, u64 sp1) {
  u32 *small = (u32 *)(uintptr_t)smallOff;
  for (u64 v = b; v >= a;) {
    u64 q = fdiv(v, p);
    u64 w = q * p;
    if (w < a) w = a;
    u32 sub = (u32)((u64)small[q] - sp1);
    sub_block(small + w, v - w + 1, sub);
    if (w == 0) break;
    v = w - 1;
  }
}
