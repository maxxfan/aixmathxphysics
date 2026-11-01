// qdfs.c -- generic lift search for an Omega-invariant F_2-linear stabilizer state (2 qubits/party, nibble per party)
// where Omega is a symplectic involution, K = S^Omega = N S (N = 1+Omega) has dim n, S = K + span{u_i}, u_i = p_i + z_i,
// N p_i = d_i (d_i = basis of K), z_i = sum_m q_im eps_m  (eps_m in V^Omega, <p_i,eps_m> = delta_im),  q_ij + q_ji = a_ij = <p_i,p_j>.
// AME <=> every word of every coset sum_i x_i p_i + sum_m c_m eps_m + K (x != 0, c = q^T x) has party weight >= THR+1.
// stdin: n kd THR  then kd hex d_i, kd hex p_i, kd hex eps_m, kd hex a_i (bit j = a_ij).  vectors are 32-bit, nibble per party.
// argv[1] = max solutions to print.  Output: SOL lines (q row masks) and "NSOL nsol nodes ...".
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int n, kd, THR;
static uint32_t d[8], p[8], eps[8], a[8];
static uint32_t Kspan[256], P[256], E[256];
static uint8_t bad[256][256];
static uint8_t row[8];
static uint8_t sums[256];
static unsigned long long nodes[9], nsol = 0, maxprint = 5, maxsol = ~0ULL;
static int stop = 0;

static inline int wt(uint32_t v) {
  uint32_t t = v | (v >> 1);
  t |= t >> 2;
  t &= 0x11111111u;
  return __builtin_popcount(t);
}

static void dfs(int k) {
  if (k == kd) {
    nsol++;
    if (nsol >= maxsol) stop = 1;
    if (nsol <= maxprint) {
      printf("SOL");
      for (int i = 0; i < kd; i++) printf(" %02x", row[i]);
      printf("\n");
      fflush(stdout);
    }
    return;
  }
  unsigned low = 0;
  for (int j = 0; j < k; j++) low |= (((row[j] >> k) & 1u) ^ ((a[j] >> k) & 1u)) << j;
  int nfree = kd - k;
  for (unsigned m = 0; m < (1u << nfree); m++) {
    unsigned c = low | (m << k);
    row[k] = (uint8_t)c;
    int ok = 1;
    for (unsigned y = 0; y < (1u << k); y++) {
      unsigned x = y | (1u << k);
      uint8_t s = sums[y] ^ (uint8_t)c;
      sums[x] = s;
      if (bad[x][s]) { ok = 0; break; }
    }
    if (!ok) continue;
    nodes[k + 1]++;
    dfs(k + 1);
    if (stop) return;
  }
}

int main(int argc, char **argv) {
  if (argc > 1) maxprint = strtoull(argv[1], 0, 10);
  if (argc > 2) maxsol = strtoull(argv[2], 0, 10);
  if (scanf("%d %d %d", &n, &kd, &THR) != 3) return 1;
  unsigned v;
  for (int i = 0; i < kd; i++) { if (scanf("%x", &v) != 1) return 1; d[i] = v; }
  for (int i = 0; i < kd; i++) { if (scanf("%x", &v) != 1) return 1; p[i] = v; }
  for (int i = 0; i < kd; i++) { if (scanf("%x", &v) != 1) return 1; eps[i] = v; }
  for (int i = 0; i < kd; i++) { if (scanf("%x", &v) != 1) return 1; a[i] = v; }
  int NX = 1 << kd;
  Kspan[0] = P[0] = E[0] = 0;
  for (int x = 1; x < NX; x++) {
    int b = __builtin_ctz(x);
    Kspan[x] = Kspan[x & (x - 1)] ^ d[b];
    P[x] = P[x & (x - 1)] ^ p[b];
    E[x] = E[x & (x - 1)] ^ eps[b];
  }
  memset(bad, 0, sizeof(bad));
  for (int x = 1; x < NX; x++)
    for (int c = 0; c < NX; c++) {
      uint32_t base = P[x] ^ E[c];
      int b = 0;
      for (int k = 0; k < NX; k++) if (wt(base ^ Kspan[k]) <= THR) { b = 1; break; }
      bad[x][c] = (uint8_t)b;
    }
  sums[0] = 0;
  dfs(0);
  printf("NSOL %llu nodes", nsol);
  for (int k = 1; k <= kd; k++) printf(" %llu", nodes[k]);
  printf("\n");
  return 0;
}
