/* Symmetric circulant graph codes over GF(q): Gamma[i][j] = c[min(|i-j|, n-|i-j|)], c[1..h], h = n/2.
   AME <=> for every A with |A| = n/2 (floor), the |A| x (n-|A|) block Gamma[A, A^c] has rank |A|.
   usage: circ TABLEFILE n mode a b seed
     TABLEFILE: q, then add table (q*q ints), then mul table (q*q ints)
     mode 0: exhaustive over c-vectors with index in [a,b) (index = base-q number, c[1] least significant), first nonzero = 1
     mode 1: random: b samples, seed
   prints "HIT c1 c2 ..." for each hit; at end "DONE tested hits".
   Move-to-front on failing subsets. Subsets reduced modulo the dihedral group. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int q, n, h, k;
static int ADD[256][256], MUL[256][256], INV[256], NEG[256];
static int nsub;
static uint32_t *subs;     /* canonical subsets */
static int (*Aidx)[24];    /* A list */
static int (*Bidx)[24];    /* complement list */
static int G[24][24];

static uint32_t rot(uint32_t m, int r) {
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    return ((m << r) | (m >> (n - r))) & full;
}
static uint32_t refl(uint32_t m) {
    uint32_t o = 0;
    for (int i = 0; i < n; i++) if (m >> i & 1) o |= 1u << ((n - i) % n);
    return o;
}

static void build_subsets(void) {
    int cap = 1 << 20;
    subs = malloc(cap * sizeof(uint32_t));
    nsub = 0;
    for (uint32_t m = 0; m < (1u << n); m++) {
        if (__builtin_popcount(m) != k) continue;
        int canon = 1;
        for (int r = 0; r < n && canon; r++) {
            uint32_t a = rot(m, r), b = refl(a);
            if (a < m || b < m) canon = 0;
        }
        if (canon) subs[nsub++] = m;
    }
    Aidx = malloc(nsub * sizeof(*Aidx));
    Bidx = malloc(nsub * sizeof(*Bidx));
    for (int s = 0; s < nsub; s++) {
        int na = 0, nb = 0;
        for (int i = 0; i < n; i++) {
            if (subs[s] >> i & 1) Aidx[s][na++] = i; else Bidx[s][nb++] = i;
        }
    }
}

/* rank of k x (n-k) block == k ? */
static int fullrank(int s) {
    int M[24][24];
    int nc = n - k;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < nc; j++) M[i][j] = G[Aidx[s][i]][Bidx[s][j]];
    int r = 0;
    for (int c = 0; c < nc && r < k; c++) {
        int pr = -1;
        for (int i = r; i < k; i++) if (M[i][c]) { pr = i; break; }
        if (pr < 0) continue;
        if (pr != r) for (int j = c; j < nc; j++) { int t = M[r][j]; M[r][j] = M[pr][j]; M[pr][j] = t; }
        int iv = INV[M[r][c]];
        for (int i = r + 1; i < k; i++) {
            if (!M[i][c]) continue;
            int f = MUL[M[i][c]][iv];
            for (int j = c; j < nc; j++)
                M[i][j] = ADD[M[i][j]][NEG[MUL[f][M[r][j]]]];
        }
        r++;
    }
    return r == k;
}

static int test(const int *c, int *firstfail) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int d = abs(i - j); if (d > n - d) d = n - d;
            G[i][j] = d ? c[d] : 0;
        }
    for (int s = 0; s < nsub; s++) {
        if (!fullrank(s)) {
            if (s > 0) { /* move to front (swap) */
                uint32_t t = subs[s]; subs[s] = subs[0]; subs[0] = t;
                int ta[24]; memcpy(ta, Aidx[s], sizeof ta); memcpy(Aidx[s], Aidx[0], sizeof ta); memcpy(Aidx[0], ta, sizeof ta);
                memcpy(ta, Bidx[s], sizeof ta); memcpy(Bidx[s], Bidx[0], sizeof ta); memcpy(Bidx[0], ta, sizeof ta);
            }
            if (firstfail) *firstfail = s;
            return 0;
        }
    }
    return 1;
}

static uint64_t rng_s[2];
static uint64_t rnd(void) {
    uint64_t s1 = rng_s[0], s0 = rng_s[1];
    rng_s[0] = s0; s1 ^= s1 << 23;
    rng_s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return rng_s[1] + s0;
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage\n"); return 1; }
    FILE *f = fopen(argv[1], "r");
    if (fscanf(f, "%d", &q) != 1) return 1;
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) if (fscanf(f, "%d", &ADD[a][b]) != 1) return 1;
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) if (fscanf(f, "%d", &MUL[a][b]) != 1) return 1;
    fclose(f);
    for (int a = 0; a < q; a++) {
        for (int b = 0; b < q; b++) { if (MUL[a][b] == 1) INV[a] = b; if (ADD[a][b] == 0) NEG[a] = b; }
    }
    n = atoi(argv[2]); int mode = atoi(argv[3]);
    long long a = atoll(argv[4]), b = atoll(argv[5]);
    uint64_t seed = strtoull(argv[6], 0, 10);
    rng_s[0] = seed * 0x9E3779B97F4A7C15ull + 1; rng_s[1] = seed ^ 0xD1B54A32D192ED03ull;
    for (int i = 0; i < 20; i++) rnd();
    h = n / 2; k = n / 2;
    build_subsets();
    fprintf(stderr, "q=%d n=%d h=%d k=%d subsets=%d\n", q, n, h, k, nsub);
    long long tested = 0, hits = 0;
    int c[24];
    if (mode == 0) {
        for (long long idx = a; idx < b; idx++) {
            long long t = idx; int first = 0;
            for (int j = 1; j <= h; j++) { c[j] = t % q; t /= q; }
            for (int j = 1; j <= h; j++) if (c[j]) { first = c[j]; break; }
            if (first != 1) continue;
            tested++;
            if (test(c, 0)) {
                hits++;
                printf("HIT"); for (int j = 1; j <= h; j++) printf(" %d", c[j]); printf("\n"); fflush(stdout);
            }
        }
    } else {
        for (long long it = 0; it < b; it++) {
            for (int j = 1; j <= h; j++) c[j] = rnd() % q;
            tested++;
            if (test(c, 0)) {
                hits++;
                printf("HIT"); for (int j = 1; j <= h; j++) printf(" %d", c[j]); printf("\n"); fflush(stdout);
            }
        }
    }
    printf("DONE tested=%lld hits=%lld\n", tested, hits);
    return 0;
}
