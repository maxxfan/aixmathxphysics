/* Auditor phase 2: independent count for the V-decomposition picture of stabilizer AME(6,4).
   stdin: nreps, then per rep 8 ints = basis r[0..7] of R = V_0^perp (12-bit vectors, party p = bits 4p..4p+3).
   stdout: per rep t = #{4-dim V_1 subset R : V_1 nondegenerate, V_2 = R cap V_1^perp, and V_1, V_2 each project
   isomorphically onto every party W_p (p=0,1,2)}.  (V_j cap (W_l+W_l') = 0 <=> projection to the third party is injective.)
   cc -O2 -o vcount vcount.c */
#include <stdio.h>
#include <stdlib.h>

static inline int sw(int v) { return ((v >> 2) & 0x333) | ((v & 0x333) << 2); }
static inline int om(int u, int v) { return __builtin_popcount(u & sw(v)) & 1; }

static int nonsing4(int a0, int a1, int a2, int a3) {
    int r[4] = {a0, a1, a2, a3};
    int rank = 0;
    for (int bit = 0; bit < 4; bit++) {
        int piv = -1;
        for (int i = rank; i < 4; i++) if ((r[i] >> bit) & 1) { piv = i; break; }
        if (piv < 0) return 0;
        int t = r[piv]; r[piv] = r[rank]; r[rank] = t;
        for (int i = 0; i < 4; i++) if (i != rank && ((r[i] >> bit) & 1)) r[i] ^= r[rank];
        rank++;
    }
    return 1;
}

static int proj_ok(const int *v) {
    for (int p = 0; p < 3; p++) {
        if (!nonsing4((v[0] >> (4 * p)) & 15, (v[1] >> (4 * p)) & 15, (v[2] >> (4 * p)) & 15, (v[3] >> (4 * p)) & 15))
            return 0;
    }
    return 1;
}

int main(void) {
    int nreps;
    if (scanf("%d", &nreps) != 1) return 1;
    for (int rep = 0; rep < nreps; rep++) {
        int r[8];
        for (int i = 0; i < 8; i++) if (scanf("%d", &r[i]) != 1) return 1;
        int omrow[8];
        for (int a = 0; a < 8; a++) {
            omrow[a] = 0;
            for (int b = 0; b < 8; b++) if (om(r[a], r[b])) omrow[a] |= 1 << b;
        }
        long long t = 0, total = 0;
        for (int mask = 0; mask < 256; mask++) {
            if (__builtin_popcount(mask) != 4) continue;
            int c[4], k = 0;
            for (int j = 0; j < 8; j++) if ((mask >> j) & 1) c[k++] = j;
            int frow[32], fcol[32], f = 0;
            for (int i = 0; i < 4; i++)
                for (int j = c[i] + 1; j < 8; j++)
                    if (!((mask >> j) & 1)) { frow[f] = i; fcol[f] = j; f++; }
            for (long long cnt = 0; cnt < (1LL << f); cnt++) {
                int row[4];
                for (int i = 0; i < 4; i++) row[i] = 1 << c[i];
                for (int b = 0; b < f; b++) if ((cnt >> b) & 1) row[frow[b]] |= 1 << fcol[b];
                total++;
                int v[4];
                for (int i = 0; i < 4; i++) {
                    int x = 0;
                    for (int a = 0; a < 8; a++) if ((row[i] >> a) & 1) x ^= r[a];
                    v[i] = x;
                }
                int g[4];
                for (int i = 0; i < 4; i++) {
                    g[i] = 0;
                    for (int j = 0; j < 4; j++) if (om(v[i], v[j])) g[i] |= 1 << j;
                }
                if (!nonsing4(g[0], g[1], g[2], g[3])) continue;
                if (!proj_ok(v)) continue;
                /* V_2 coordinates: null space of rows M_i = sum_{a in row_i} omrow[a] */
                int M[4];
                for (int i = 0; i < 4; i++) {
                    int x = 0;
                    for (int a = 0; a < 8; a++) if ((row[i] >> a) & 1) x ^= omrow[a];
                    M[i] = x;
                }
                int piv[4], rk = 0;
                for (int col = 0; col < 8 && rk < 4; col++) {
                    int p = -1;
                    for (int i = rk; i < 4; i++) if ((M[i] >> col) & 1) { p = i; break; }
                    if (p < 0) continue;
                    int tmp = M[p]; M[p] = M[rk]; M[rk] = tmp;
                    for (int i = 0; i < 4; i++) if (i != rk && ((M[i] >> col) & 1)) M[i] ^= M[rk];
                    piv[rk++] = col;
                }
                if (rk != 4) { fprintf(stderr, "rank deficiency\n"); return 2; }
                int ispiv = 0;
                for (int i = 0; i < 4; i++) ispiv |= 1 << piv[i];
                int w[4], wn = 0;
                for (int col = 0; col < 8; col++) {
                    if ((ispiv >> col) & 1) continue;
                    int y = 1 << col;
                    for (int i = 0; i < 4; i++) if ((M[i] >> col) & 1) y |= 1 << piv[i];
                    int x = 0;
                    for (int a = 0; a < 8; a++) if ((y >> a) & 1) x ^= r[a];
                    w[wn++] = x;
                }
                if (!proj_ok(w)) continue;
                t++;
            }
        }
        printf("%lld %lld\n", t, total);
        fflush(stdout);
    }
    return 0;
}
