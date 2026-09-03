/* engine_par.c — range kernels for the multi-core engine (parallel.js).
 * Line-for-line the loops of engine.c, restricted to an index range so K
 * threads can each run their slice over ONE shared linear memory.  Built
 * with atomics + shared memory (see wasmbuild.js); the barriers live in JS.
 *
 * The `large` table is stored as IEEE doubles (exact below 2^53, and every
 * value here is ≤ x < 9×10^15) so that the JavaScript coordinator can read
 * and write the same bytes through a Float64Array; `small` is u32 on both
 * sides.  Mixing u64 and f64 views of one buffer was a bug we caught.
 */
#include <stdint.h>

typedef uint64_t u64;
typedef uint32_t u32;

/* exact floor(a/b) for a < 2^53 via pipelined double division (see engine.c) */
static inline u64 fdiv(u64 a, u64 b) {
  return (u64)((double)a / (double)b);
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

__attribute__((export_name("large_range")))
void large_range(u32 smallOff, u32 largeOff, u64 a, u64 b, u64 p, u64 sp1, u64 xp, u64 iSw) {
  u32 *small = (u32 *)(uintptr_t)smallOff;
  double *large = (double *)(uintptr_t)largeOff;
  double s = (double)sp1;
  u64 e1 = b < iSw ? b : iSw;
  u64 i = a;
  for (; i <= e1; i++) large[i] -= large[i * p] - s;
  for (; i <= b; i++) large[i] -= (double)small[fdiv(xp, i)] - s;
}

__attribute__((export_name("small_range")))
void small_range(u32 smallOff, u64 a, u64 b, u64 p, u64 sp1) {
  u32 *small = (u32 *)(uintptr_t)smallOff;
  for (u64 v = b; v >= a;) {
    u64 q = fdiv(v, p);
    u64 w = q * p;
    if (w < a) w = a;
    u32 sub = (u32)((u64)small[q] - sp1);
    for (; v >= w; v--) small[v] -= sub;
  }
}
