/* rankcheck: AME test for a stabilizer state given by an F_p generator matrix.
 *
 * stdin:  p N n m t
 *         N rows of 2N integers in [0,p)   (row = X part on N qudits, then Z part)
 * Qudit i belongs to party i / m (m qudits of dimension p per party, local dimension p^m).
 * For every t-subset B of the n parties (lexicographic order, as itertools.combinations), the
 * N x 2mt submatrix on the X and Z columns of B's qudits must have rank N over F_p; a deficient B
 * means some nonzero stabilizer is supported on the complement of B.
 *
 * stdout: "TOTAL <C(n,t)> BAD <count>"
 *         "MAP <hex>"  bit i (LSB-first within each hex digit, digits in order) = subset i deficient
 * Exit 2 on malformed input. Pure integer arithmetic; p prime <= 251 (checked by the caller).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 192   /* rows (qudits) */
#define MAXC 384   /* columns of a submatrix: 2*m*t <= 2*N */
#define MAXW 6     /* 64-bit words for p = 2 rows */

static int p, N, n, m, t;
static uint8_t S[MAXN][2 * MAXN];
static int inv[256];

static int rank_p2(const int *cols, int nc) {
    uint64_t R[MAXN][MAXW];
    int W = (nc + 63) / 64;
    for (int i = 0; i < N; i++) {
        memset(R[i], 0, sizeof(R[i]));
        for (int k = 0; k < nc; k++)
            if (S[i][cols[k]]) R[i][k >> 6] |= (uint64_t)1 << (k & 63);
    }
    int rank = 0;
    for (int c = 0; c < nc && rank < N; c++) {
        if (nc - c < N - rank) break;
        int w = c >> 6; uint64_t b = (uint64_t)1 << (c & 63);
        int piv = -1;
        for (int i = rank; i < N; i++) if (R[i][w] & b) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != rank) for (int k = 0; k < W; k++) { uint64_t x = R[piv][k]; R[piv][k] = R[rank][k]; R[rank][k] = x; }
        for (int i = rank + 1; i < N; i++)
            if (R[i][w] & b) for (int k = w; k < W; k++) R[i][k] ^= R[rank][k];
        rank++;
    }
    return rank;
}

static int rank_p(const int *cols, int nc) {
    static uint16_t R[MAXN][MAXC];
    for (int i = 0; i < N; i++)
        for (int k = 0; k < nc; k++) R[i][k] = S[i][cols[k]];
    int rank = 0;
    for (int c = 0; c < nc && rank < N; c++) {
        if (nc - c < N - rank) break;
        int piv = -1;
        for (int i = rank; i < N; i++) if (R[i][c]) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != rank) for (int k = c; k < nc; k++) { uint16_t x = R[piv][k]; R[piv][k] = R[rank][k]; R[rank][k] = x; }
        int iv = inv[R[rank][c]];
        for (int k = c; k < nc; k++) R[rank][k] = (uint16_t)((R[rank][k] * iv) % p);
        for (int i = rank + 1; i < N; i++) {
            int f = R[i][c];
            if (!f) continue;
            for (int k = c; k < nc; k++) R[i][k] = (uint16_t)((R[i][k] + (p - f) * R[rank][k]) % p);
        }
        rank++;
    }
    return rank;
}

int main(void) {
    if (scanf("%d %d %d %d %d", &p, &N, &n, &m, &t) != 5) return 2;
    if (p < 2 || p > 251 || N < 1 || N > MAXN || n < 2 || m < 1 || N != n * m || t < 1 || t > n) return 2;
    if (2 * m * t > MAXC || (p == 2 && 2 * m * t > 64 * MAXW)) return 2;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 2 * N; j++) {
            int x;
            if (scanf("%d", &x) != 1 || x < 0 || x >= p) return 2;
            S[i][j] = (uint8_t)x;
        }
    for (int a = 1; a < p; a++)
        for (int b = 1; b < p; b++)
            if ((a * b) % p == 1) inv[a] = b;
    int idx[64], cols[MAXC];
    for (int i = 0; i < t; i++) idx[i] = i;
    size_t cap = 1 << 20, total = 0, bad = 0;
    uint8_t *map = calloc(cap, 1);
    if (!map) return 3;
    for (;;) {
        int nc = 0;
        for (int i = 0; i < t; i++) {
            int j = idx[i];
            for (int k = 0; k < m; k++) cols[nc++] = j * m + k;
            for (int k = 0; k < m; k++) cols[nc++] = N + j * m + k;
        }
        int r = (p == 2) ? rank_p2(cols, nc) : rank_p(cols, nc);
        if (total / 8 >= cap) {
            map = realloc(map, cap * 2);
            if (!map) return 3;
            memset(map + cap, 0, cap);
            cap *= 2;
        }
        if (r < N) { bad++; map[total / 8] |= (uint8_t)(1 << (total % 8)); }
        total++;
        int i = t - 1;
        while (i >= 0 && idx[i] == n - t + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int k = i + 1; k < t; k++) idx[k] = idx[k - 1] + 1;
    }
    printf("TOTAL %zu BAD %zu\nMAP ", total, bad);
    static const char hx[] = "0123456789abcdef";
    for (size_t b = 0; b < (total + 3) / 4; b++) {
        int nib = (map[b / 2] >> ((b % 2) * 4)) & 15;
        putchar(hx[nib]);
    }
    putchar('\n');
    free(map);
    return 0;
}
