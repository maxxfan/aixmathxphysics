/* cutrank.c -- release checker, written independently of the project verifier.
 *
 * Input (stdin): p N n m t, then N rows of 2N integers in [0,p) (X part then Z part; qudit i belongs to party i/m).
 * For every t-subset B of {0..n-1} (lexicographic order) it computes the rank over F_p of the N x 2mt column block of
 * the X and Z columns of the qudits of B, and counts the subsets whose rank is < N.
 * Output: "SETS <C(n,t)> DEFICIENT <count>" and, if any, up to 20 lines "DEF i0 i1 ...".
 * Rank: Gaussian elimination by columns; p = 2 uses 64-bit words (up to 4 words), odd p uses int arithmetic mod p.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXR 128
#define MAXC 256

static int P, NR, NP, M, T;
static int A[MAXR][2 * MAXR];

static int inv_mod(int a) { for (int x = 1; x < P; x++) if ((a * x) % P == 1) return x; return 0; }

static int rank_block(const int *cols, int nc) {
    if (P == 2) {
        uint64_t row[MAXR][4];
        int W = (nc + 63) / 64;
        for (int r = 0; r < NR; r++) {
            for (int w = 0; w < 4; w++) row[r][w] = 0;
            for (int c = 0; c < nc; c++) if (A[r][cols[c]]) row[r][c >> 6] |= (uint64_t)1 << (c & 63);
        }
        int rk = 0;
        for (int c = 0; c < nc && rk < NR; c++) {
            int w = c >> 6; uint64_t bit = (uint64_t)1 << (c & 63);
            int pv = -1;
            for (int r = rk; r < NR; r++) if (row[r][w] & bit) { pv = r; break; }
            if (pv < 0) continue;
            for (int x = 0; x < W; x++) { uint64_t t = row[pv][x]; row[pv][x] = row[rk][x]; row[rk][x] = t; }
            for (int r = 0; r < NR; r++) if (r != rk && (row[r][w] & bit)) for (int x = 0; x < W; x++) row[r][x] ^= row[rk][x];
            rk++;
        }
        return rk;
    } else {
        static int B[MAXR][MAXC];
        for (int r = 0; r < NR; r++) for (int c = 0; c < nc; c++) B[r][c] = A[r][cols[c]];
        int rk = 0;
        for (int c = 0; c < nc && rk < NR; c++) {
            int pv = -1;
            for (int r = rk; r < NR; r++) if (B[r][c]) { pv = r; break; }
            if (pv < 0) continue;
            for (int x = 0; x < nc; x++) { int t = B[pv][x]; B[pv][x] = B[rk][x]; B[rk][x] = t; }
            int iv = inv_mod(B[rk][c]);
            for (int x = 0; x < nc; x++) B[rk][x] = (B[rk][x] * iv) % P;
            for (int r = 0; r < NR; r++) {
                if (r == rk || !B[r][c]) continue;
                int f = B[r][c];
                for (int x = 0; x < nc; x++) B[r][x] = ((B[r][x] - f * B[rk][x]) % P + P) % P;
            }
            rk++;
        }
        return rk;
    }
}

int main(void) {
    if (scanf("%d %d %d %d %d", &P, &NR, &NP, &M, &T) != 5) return 2;
    if (P < 2 || P > 251 || NR != NP * M || NR > MAXR || T < 1 || T > NP || 2 * M * T > MAXC) return 2;
    if (P == 2 && 2 * M * T > 256) return 2;
    for (int i = 0; i < NR; i++) for (int j = 0; j < 2 * NR; j++) {
        if (scanf("%d", &A[i][j]) != 1 || A[i][j] < 0 || A[i][j] >= P) return 2;
    }
    int idx[64], cols[MAXC];
    for (int i = 0; i < T; i++) idx[i] = i;
    long long sets = 0, def = 0;
    for (;;) {
        int nc = 0;
        for (int i = 0; i < T; i++) {
            for (int k = 0; k < M; k++) cols[nc++] = idx[i] * M + k;
            for (int k = 0; k < M; k++) cols[nc++] = NR + idx[i] * M + k;
        }
        if (rank_block(cols, nc) < NR) {
            if (def < 20) { printf("DEF"); for (int i = 0; i < T; i++) printf(" %d", idx[i]); printf("\n"); }
            def++;
        }
        sets++;
        int i = T - 1;
        while (i >= 0 && idx[i] == NP - T + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i + 1; j < T; j++) idx[j] = idx[j - 1] + 1;
    }
    printf("SETS %lld DEFICIENT %lld\n", sets, def);
    return 0;
}
