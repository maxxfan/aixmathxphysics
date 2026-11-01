// inv8.c -- AME(8,4) stabilizer states (2 qubits/party) invariant under a fixed-point-free party involution
// (party p <-> p+4, trivial local lift, WLOG).  Structure (pairs shift 3, see NOTES):
//   S^sigma = N S = C (N = 1+sigma) folds to a Lagrangian C in W = (F_2^4)^4 with all nonzero words of position-weight >= 3
//   (an AME(4,4) stabilizer state), and S = span{ (u_i, u_i + d_i), (d_i, d_i) : i<8 }, d_i = basis of C,
//   u_i = sum_j q_ij e_j (e_j = pure-Z dual basis), q symmetric 8x8 over F_2 (36 free bits).
//   Word with coefficient x (nonzero): (U, U+D(x)), U in coset of C labelled by kappa = q x (kappa_k = <U, d_k>).
//   AME  <=>  for all x != 0:  q x not in B_x  (B_x = cosets containing U with wt(U)+wt(U+D) <= 4).
// W vector: bits 0..7 = x part (qubit b), bits 8..15 = z part; position p = qubits 2p,2p+1.
// Input (stdin): 8 hex ints d_0..d_7.   argv[1] = max solutions to print (default 5).   argv[2] = (optional) max total
// solutions to count before stopping (default unlimited).
// Output: "NSOL count nodes per level ..." and "SOL c0 .. c7" (q columns as 8-bit masks, bit j of c_i = q_ij).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static uint8_t wtW[65536];
static uint16_t d[8];
static uint8_t bad[256][256];
static uint8_t col[8];
static uint8_t sums[256];
static unsigned long long nodes[9], nsol = 0, maxprint = 5, maxsol = ~0ULL;
static int stop = 0, THR = 4;  // word is bad if folded weight wt(U)+wt(U+D) <= THR (4 = AME(8,4); 3 = distance-4 control)

static inline int sp(unsigned a, unsigned b) {
  return __builtin_parity(((a & 255u) & (b >> 8)) ^ ((a >> 8) & (b & 255u)));
}

static int DUMP7 = 0;

static void dfs(int k) {
  if (stop) return;
  if (k == 7 && DUMP7) {
    printf("L7");
    for (int i = 0; i < 7; i++) printf(" %02x", col[i]);
    printf("\n");
  }
  if (k == 8) {
    nsol++;
    if (nsol <= maxprint) {
      printf("SOL");
      for (int i = 0; i < 8; i++) printf(" %02x", col[i]);
      printf("\n");
      fflush(stdout);
    }
    if (nsol >= maxsol) stop = 1;
    return;
  }
  int nfree = 8 - k;
  unsigned lowmask = 0;
  for (int j = 0; j < k; j++) lowmask |= ((col[j] >> k) & 1u) << j;  // q_{j,k} = q_{k,j}
  for (unsigned m = 0; m < (1u << nfree); m++) {
    unsigned c = lowmask | (m << k);
    col[k] = (uint8_t)c;
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
  if (argc > 3) THR = atoi(argv[3]);
  if (argc > 4) DUMP7 = atoi(argv[4]);
  for (int w = 0; w < 65536; w++) {
    int c = 0;
    for (int p = 0; p < 4; p++) if (w & ((3u << (2 * p)) | (3u << (8 + 2 * p)))) c++;
    wtW[w] = (uint8_t)c;
  }
  for (int i = 0; i < 8; i++) { unsigned v; if (scanf("%x", &v) != 1) return 1; d[i] = (uint16_t)v; }
  memset(bad, 0, sizeof(bad));
  for (unsigned x = 1; x < 256; x++) {
    unsigned D = 0;
    for (int i = 0; i < 8; i++) if (x >> i & 1) D ^= d[i];
    for (unsigned U = 0; U < 65536; U++) {
      if (wtW[U] + wtW[U ^ D] <= THR) {
        unsigned kap = 0;
        for (int k = 0; k < 8; k++) kap |= (unsigned)sp(U, d[k]) << k;
        bad[x][kap] = 1;
      }
    }
  }
  sums[0] = 0;
  dfs(0);
  printf("NSOL %llu nodes", nsol);
  for (int k = 1; k <= 8; k++) printf(" %llu", nodes[k]);
  printf(" stopped=%d\n", stop);
  return 0;
}
