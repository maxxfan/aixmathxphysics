/* lagsearch2: staged search.  L (dim 6, Lagrangian in V^Omega = (a,b,c,d), codes: bits 0-3 a, 4-7 b, 8-9 c, 10-11 d) must have
 * L_a0 = L cap {a=0} of dim exactly 2 (words need b,c,d != 0), L_b0 = L cap {b=0} of dim exactly 2, L_a0 cap L_b0 = 0.
 * So basis: v1,v2 in Va0, v3,v4 in Vb0, v5,v6 arbitrary in Vs.  usage: lagsearch2 FILE FIRST LAST */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static unsigned char inV[4096];
static int S[3][4096], nS[3];
static long nodes, nfound;
static int firstL[6], haveL;

static inline int om(int u, int v) {
    int sv = ((v & 0x555) << 1) | ((v >> 1) & 0x555);
    return __builtin_popcount(u & sv) & 1;
}

static void rec(int k, int *basis, int *U, int nU) {
    nodes++;
    if (k == 6) {
        nfound++;
        if (!haveL) { memcpy(firstL, basis, sizeof(int) * 6); haveL = 1; }
        return;
    }
    int st = k / 2;
    int U2[64];
    for (int i = 0; i < nS[st]; i++) {
        int v = S[st][i];
        if ((k & 1) && v <= basis[k - 1]) continue;
        int ok = 1;
        for (int j = 0; j < nU; j++) if (U[j] == v) { ok = 0; break; }
        if (!ok) continue;
        for (int j = 0; j < k; j++) if (om(v, basis[j])) { ok = 0; break; }
        if (!ok) continue;
        for (int j = 0; j < nU; j++) if (!inV[v ^ U[j]]) { ok = 0; break; }
        if (!ok) continue;
        for (int j = 0; j < nU; j++) { U2[j] = U[j]; U2[nU + j] = U[j] ^ v; }
        basis[k] = v;
        rec(k + 1, basis, U2, nU * 2);
    }
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb");
    int first = atoi(argv[2]), last = atoi(argv[3]);
    int nT; fread(&nT, 4, 1, f);
    for (int t = 0; t < nT; t++) {
        uint16_t N; fread(&N, 2, 1, f);
        uint16_t *cs = (uint16_t *)malloc(2 * N); fread(cs, 2, N, f);
        if (t < first || t > last) { free(cs); continue; }
        memset(inV, 0, sizeof inV);
        nS[0] = nS[1] = nS[2] = 0;
        for (int i = 0; i < N; i++) {
            int c = cs[i]; inV[c] = 1;
            int a = c & 15, b = (c >> 4) & 15;
            if (c == 0) continue;
            if (a == 0) S[0][nS[0]++] = c;
            if (b == 0) S[1][nS[1]++] = c;
            S[2][nS[2]++] = c;
        }
        nodes = 0; nfound = 0; haveL = 0;
        int basis[8]; int U[1] = {0};
        rec(0, basis, U, 1);
        printf("T %d |Vs| %d |Va0| %d |Vb0| %d nodes %ld Lagrangians %ld", t, N, nS[0], nS[1], nodes, nfound);
        if (haveL) { printf(" first:"); for (int i = 0; i < 6; i++) printf(" %d", firstL[i]); }
        printf("\n"); fflush(stdout);
        free(cs);
    }
    return 0;
}
