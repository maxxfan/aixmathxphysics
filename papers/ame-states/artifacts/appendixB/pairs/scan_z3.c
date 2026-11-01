// scan_z3.c : sector sigma = (0 1 2)(3 4 5), fixed parties 6,7.  S = (lift of an AME(4,4) state C_1 on positions (A,B,f1,f2)) + S_nf,
// S_nf = graph of an F_4 anti-isometry V_A -> V_B.  Input file (argv[1]): "R M" then R lines of 8 hex ints (C_1 lifts, 32-bit V vectors
// x | z<<16) then M lines of 8 hex ints (S_nf generators).  For every pair (r,m): test all 2^16 words, AME iff no nonzero word of party weight <= 4.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static uint32_t pm[8];

static inline int wt(uint32_t w) {
  int c = 0;
  for (int j = 0; j < 8; j++) c += (w & pm[j]) != 0;
  return c;
}

int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r");
  int R, M;
  if (fscanf(f, "%d %d", &R, &M) != 2) return 1;
  for (int j = 0; j < 8; j++) pm[j] = (3u << (2 * j)) | (3u << (16 + 2 * j));
  uint32_t (*reps)[8] = malloc(sizeof(uint32_t) * 8 * R);
  uint32_t (*phis)[8] = malloc(sizeof(uint32_t) * 8 * M);
  for (int r = 0; r < R; r++) for (int i = 0; i < 8; i++) { unsigned v; if (fscanf(f, "%x", &v) != 1) return 2; reps[r][i] = v; }
  for (int m = 0; m < M; m++) for (int i = 0; i < 8; i++) { unsigned v; if (fscanf(f, "%x", &v) != 1) return 3; phis[m][i] = v; }
  fclose(f);
  // S_nf alone
  char *okphi = calloc(M, 1);
  int nok = 0;
  for (int m = 0; m < M; m++) {
    uint32_t w = 0; int good = 1;
    for (unsigned i = 1; i < 256; i++) { w ^= phis[m][__builtin_ctz(i)]; if (wt(w) <= 4) { good = 0; break; } }
    okphi[m] = good; nok += good;
  }
  printf("R=%d M=%d phi passing S_nf alone: %d\n", R, M, nok); fflush(stdout);
  int *list = malloc(sizeof(int) * (nok + 1)); int nl = 0;
  for (int m = 0; m < M; m++) if (okphi[m]) list[nl++] = m;
  unsigned long long nsol = 0, tested = 0;
  int best = 99; // minimum number of low-weight words seen is not tracked; only AME solutions
  for (int r = 0; r < R; r++) {
    for (int li = 0; li < nl; li++) {
      int m = list[li];
      uint32_t g[16];
      for (int i = 0; i < 8; i++) { g[i] = reps[r][i]; g[8 + i] = phis[m][i]; }
      uint32_t w = 0; int good = 1;
      for (unsigned i = 1; i < 65536; i++) { w ^= g[__builtin_ctz(i)]; if (wt(w) <= 4) { good = 0; break; } }
      tested++;
      if (good) {
        nsol++;
        if (nsol <= 5) { printf("SOL r=%d m=%d\n", r, m); fflush(stdout); }
      }
    }
    if (r % 50 == 0) { printf("progress r=%d/%d tested=%llu nsol=%llu\n", r, R, tested, nsol); fflush(stdout); }
  }
  printf("DONE R=%d tested=%llu NSOL=%llu\n", R, tested, nsol);
  (void)best;
  return 0;
}
