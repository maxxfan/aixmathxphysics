// dfs4.c -- exhaustive DFS for F_4 "graph" AME(n,4) (n even or odd):
// matrix h over F_4, zero diagonal, Hermitian (type 0: h_ji = conj h_ij) or symmetric (type 1: h_ji = h_ij).
// AME  <=>  every S, |S| <= floor(n/2), has rank_F4 h[S, S^c] = |S|.
// Symmetry breaking: scaling D h D^* normalises row 0 to {0,1} with m ones (vertices 1..m), sorted; row 1 sorted
// inside classes {2..m}, {m+1..n-1} and normalised to {0,1} on the class {m+1..n-1}.
// usage: dfs4 n type m part nparts tlimit    (reads nothing)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

static int n, type, K;
static uint8_t mul[4][4], inv_[4], cj[4];
static uint8_t h[12][12];
static long nodes[12];
static long nsol = 0;
static int nsets[12];
static int sets[12][600][2]; // [level][i] = {mask, size}
static clock_t t0;
static double tlimit;
static int timedout = 0;

static void init_f4(void) {
  // elements 0,1,2=w,3=w^2=w+1 ; add = xor
  for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) {
    if (a == 0 || b == 0) mul[a][b] = 0;
    else if (a == 1) mul[a][b] = b;
    else if (b == 1) mul[a][b] = a;
    else if (a == 2 && b == 2) mul[a][b] = 3;
    else if (a == 3 && b == 3) mul[a][b] = 2;
    else mul[a][b] = 1;
  }
  inv_[1] = 1; inv_[2] = 3; inv_[3] = 2;
  cj[0] = 0; cj[1] = 1; cj[2] = 3; cj[3] = 2;
}

static int rank_set(int mask) {
  int rows[12], cols[12], r = 0, c = 0;
  for (int i = 0; i < n; i++) { if (mask >> i & 1) rows[r++] = i; else cols[c++] = i; }
  uint8_t M[6][12];
  for (int a = 0; a < r; a++) for (int b = 0; b < c; b++) M[a][b] = h[rows[a]][cols[b]];
  int rk = 0;
  for (int col = 0; col < c && rk < r; col++) {
    int pr = -1;
    for (int a = rk; a < r; a++) if (M[a][col]) { pr = a; break; }
    if (pr < 0) continue;
    if (pr != rk) for (int b = 0; b < c; b++) { uint8_t t = M[pr][b]; M[pr][b] = M[rk][b]; M[rk][b] = t; }
    uint8_t iv = inv_[M[rk][col]];
    for (int b = col; b < c; b++) M[rk][b] = mul[iv][M[rk][b]];
    for (int a = rk + 1; a < r; a++) if (M[a][col]) {
      uint8_t f = M[a][col];
      for (int b = col; b < c; b++) M[a][b] ^= mul[f][M[rk][b]];
    }
    rk++;
  }
  return rk;
}

static int check_level(int k) {
  for (int i = 0; i < nsets[k]; i++) if (rank_set(sets[k][i][0]) != sets[k][i][1]) return 0;
  return 1;
}

static void set_entry(int i, int j, int v) { h[i][j] = v; h[j][i] = type == 0 ? cj[v] : v; }

static int final_check(void) {
  for (int mask = 0; mask < (1 << n); mask++) {
    if (!(mask >> (n - 1) & 1)) continue;
    int sz = __builtin_popcount(mask);
    if (sz > K) continue;
    if (rank_set(mask) != sz) return 0;
  }
  return 1;
}

static void report(void) {
  if (!final_check()) return;
  nsol++;
  if (nsol <= 20) {
    printf("SOL");
    for (int i = 0; i < n; i++) { printf(" "); for (int j = 0; j < n; j++) printf("%d", h[i][j]); }
    printf("\n"); fflush(stdout);
  }
}

static int mm, part_, nparts_;
static long idx1 = 0;
static uint8_t same[12][12]; // same[k][j]: j and j-1 in same cell going into row k

static void dfs_row(int k, int j) {
  if (timedout) return;
  if (j == n) {
    if (k == 1) { if (idx1++ % nparts_ != part_) return; }
    nodes[k]++;
    if ((nodes[k] & 0xFFFFF) == 0 && (double)(clock() - t0) / CLOCKS_PER_SEC > tlimit) { timedout = 1; return; }
    if (!check_level(k)) return;
    if (k == n - 2) { report(); return; }
    // refine cells
    for (int jj = 0; jj < n; jj++) same[k + 1][jj] = 0;
    for (int jj = k + 2; jj < n; jj++) same[k + 1][jj] = same[k][jj] && h[k][jj - 1] == h[k][jj];
    dfs_row(k + 1, k + 2);
    return;
  }
  int vmax = (k == 1 && j > mm) ? 2 : 4;
  int lo = 0;
  if (j >= k + 2 && same[k][j]) lo = h[k][j - 1];
  for (int v = lo; v < vmax; v++) { set_entry(k, j, v); dfs_row(k, j + 1); }
}

int main(int argc, char **argv) {
  n = atoi(argv[1]); type = atoi(argv[2]); int m = atoi(argv[3]);
  part_ = atoi(argv[4]); nparts_ = atoi(argv[5]); tlimit = atof(argv[6]);
  mm = m;
  K = n / 2;
  init_f4();
  for (int k = 0; k <= n - 2; k++) {
    int c = 0;
    for (int sz = 1; sz <= K; sz++)
      for (int mask = 0; mask < (1 << (k + 1)); mask++)
        if ((mask >> k & 1) && __builtin_popcount(mask) == sz) { sets[k][c][0] = mask; sets[k][c][1] = sz; c++; }
    nsets[k] = c;
  }
  memset(h, 0, sizeof h);
  memset(same, 0, sizeof same);
  t0 = clock();
  for (int j = 1; j < n; j++) set_entry(0, j, j <= m ? 1 : 0);
  if (!check_level(0)) { printf("row0 fails\n"); return 0; }
  nodes[0] = 1;
  for (int jj = 2; jj < n; jj++) same[1][jj] = (jj <= m && jj >= 2 && jj - 1 >= 2) || (jj >= m + 2);
  // vertex 1 is the representative of the ones-cell; cell of ones among 2..m
  dfs_row(1, 2);
  printf("DONE n=%d type=%d m=%d part=%d/%d nsol=%ld timedout=%d time=%.1fs nodes:", n, type, m, part_, nparts_, nsol, timedout, (double)(clock() - t0) / CLOCKS_PER_SEC);
  for (int k = 0; k <= n - 2; k++) printf(" %ld", nodes[k]);
  printf("\n");
  return 0;
}
