/* lagsearch: for each T read from the .bin (z3twist4.py), count Lagrangians (isotropic subspaces of dim 6) of F_2^12 (omega: pair-swap form)
 * all of whose nonzero words lie in Vs.  usage: lagsearch FILE FIRST LAST   (T indices)  prints per-T counts and the first Lagrangian basis.
 * Canonical enumeration: basis v1<v2<..., each v_i the minimum of its coset mod span(previous). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int ldim = 6;
static unsigned char inV[4096];
static long nodes, nfound;
static int firstL[6], haveL;

static inline int om(int u, int v) {
    int sv = ((v & 0x555) << 1) | ((v >> 1) & 0x555);
    return __builtin_popcount(u & sv) & 1;
}

static void rec(int k, int *basis, int *U, int nU, int *pool, int npool) {
    nodes++;
    if (k == ldim) {
        nfound++;
        if (!haveL) { memcpy(firstL, basis, sizeof(int) * ldim); haveL = 1; }
        return;
    }
    if (npool < ldim - k) return;
    int *npoolbuf = (int *)malloc(sizeof(int) * (npool + 1));
    int *U2 = (int *)malloc(sizeof(int) * nU * 2);
    for (int i = 0; i < npool; i++) {
        int v = pool[i];
        /* v must be the min of its coset mod U : v < v^u for all u in U */
        int ok = 1;
        for (int j = 1; j < nU; j++) if ((v ^ U[j]) < v) { ok = 0; break; }
        if (!ok) continue;
        /* new subspace */
        for (int j = 0; j < nU; j++) { U2[j] = U[j]; U2[nU + j] = U[j] ^ v; }
        int np = 0;
        for (int m = i + 1; m < npool; m++) {
            int w = pool[m];
            if (om(w, v)) continue;
            int good = 1;
            for (int j = 0; j < nU; j++) if (!inV[w ^ U[j] ^ v]) { good = 0; break; }
            if (good) npoolbuf[np++] = w;
        }
        basis[k] = v;
        /* also all words v^u must be in V (they are: v in pool is checked against U before) */
        rec(k + 1, basis, U2, nU * 2, npoolbuf, np);
    }
    free(npoolbuf); free(U2);
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
        for (int i = 0; i < N; i++) inV[cs[i]] = 1;
        int *pool = (int *)malloc(sizeof(int) * N);
        int np = 0;
        for (int i = 0; i < N; i++) pool[np++] = cs[i];
        /* sort ascending */
        for (int i = 1; i < np; i++) { int x = pool[i], j = i - 1; while (j >= 0 && pool[j] > x) { pool[j + 1] = pool[j]; j--; } pool[j + 1] = x; }
        nodes = 0; nfound = 0; haveL = 0;
        int basis[8]; int U[1] = {0};
        /* U must contain 0 and enumerates words; start with U={0}; words check uses U[j]^v ... */
        rec(0, basis, U, 1, pool, np);
        printf("T %d |Vs| %d nodes %ld Lagrangians %ld", t, N, nodes, nfound);
        if (haveL) { printf(" first:"); for (int i = 0; i < ldim; i++) printf(" %d", firstL[i]); }
        printf("\n"); fflush(stdout);
        free(cs); free(pool);
    }
    return 0;
}
