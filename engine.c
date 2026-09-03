/* engine.c — Lucy_Hedgehog exact prime counting in freestanding C,
 * compiled to WebAssembly (see wasmbuild.js).  Direct port of the verified
 * primeCount() in nthprime.js; the test suite requires digit-for-digit
 * agreement between this, the JS Lucy engine, and the JS LMO engine.
 *
 * No libc.  Memory is the module's exported linear memory; the caller
 * grows it and passes byte offsets:
 *   small: u32[r+1] at smallOff,  large: u64[r+1] at largeOff,  r = isqrt(x)
 *
 * Build: clang --target=wasm32 -O3 -nostdlib (see wasmbuild.js)
 */
#include <stdint.h>

typedef uint64_t u64;
typedef uint32_t u32;

/* small[0..len) -= sub — a run of equal updates.  (A SIMD version was
 * measured at only +1–3%: this algorithm is memory-bound, not ALU-bound,
 * so it was dropped to keep the widest browser compatibility.) */
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

/* floor(a/b) via pipelined double division — exact for a < 2^53 (the IEEE
 * proof in nthprime.js applies verbatim; every caller stays below 9e15).
 * ~5x lower latency than i64.div_u on common hardware. */
static inline u64 fdiv(u64 a, u64 b) {
  return (u64)((double)a / (double)b);
}

__attribute__((export_name("isqrt")))
u64 isqrt_export(u64 n) { return isqrt64(n); }

/* progress callback into the host (fraction = p / r), every ~2048 primes */
__attribute__((import_name("progress"))) void host_progress(u64 p, u64 r);

__attribute__((export_name("pi_lucy")))
u64 pi_lucy(u64 x, u32 smallOff, u32 largeOff) {
  if (x < 2) return 0;
  u64 r = isqrt64(x);
  u32 *small = (u32 *)(uintptr_t)smallOff;
  u64 *large = (u64 *)(uintptr_t)largeOff;

  for (u64 v = 1; v <= r; v++) small[v] = (u32)(v - 1);
  for (u64 i = 1; i <= r; i++) large[i] = fdiv(x, i) - 1;

  for (u64 p = 2; p <= r; p++) {
    if ((p & 2047) == 0) host_progress(p, r);
    if (small[p] == small[p - 1]) continue; /* p composite */
    u64 sp1 = small[p - 1];                 /* pi(p-1) */
    u64 p2 = p * p;
    u64 xp = fdiv(x, p);
    u64 imax = fdiv(x, p2);
    if (imax > r) imax = r;
    u64 isw = r / p;
    if (isw > imax) isw = imax;
    u64 i = 1;
    for (; i <= isw; i++) large[i] -= large[i * p] - sp1;
    for (; i <= imax; i++) large[i] -= (u64)small[fdiv(xp, i)] - sp1;
    if (p2 <= r) {
      /* descending blocked update: small[q] still holds the previous
       * pass's value because q = v/p < w <= v for every v in the block */
      for (u64 v = r; v >= p2;) {
        u64 q = fdiv(v, p);
        u32 sub = (u32)((u64)small[q] - sp1);
        u64 w = q * p;
        if (w < p2) w = p2;
        sub_block(small + w, v - w + 1, sub);
        v = w - 1;
      }
    }
  }
  return large[1];
}
