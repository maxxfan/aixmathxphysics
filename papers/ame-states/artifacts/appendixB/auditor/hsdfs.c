/* hsdfs.c -- auditor's independent exhaustive search for S-type AME(n,4), n = 2h parties:
 * F_4-LINEAR stabilizer with one F_4^2 = (a|b) per party (ququart as F_4 Pauli X^a Z^b), self-dual for the F_4-bilinear
 * alternating form J = sum_c (a_c b'_c + a'_c b_c).  S = rowspace [I | M], 2h x 2h, parties = coordinate pairs (2j, 2j+1);
 * self-duality  <=>  M J M^T = J  (rows of M form hyperbolic frames: <r_2p, r_2p+1> = 1, all other pairs 0).
 * AME <=> every party block minor M[R,C], |R|=|C|=k, is nonsingular <=> all C(2h,h) column sets of G=[I|M] have rank 2h.
 * Normal form: right multiplication by block-diag SL(2,4) (order 60) fixes row party 0 up to the determinants delta_c of its
 * 2x2 blocks: X_c = diag(1, delta_c), sum delta_c = 1, delta_c in F_4^*.  Later parties: canonical hyperbolic frame
 * (u = smallest nonzero vector of the plane, w = smallest vector of the plane with <u,w> = 1).
 * #solutions = N_norm * 60^h * 60^(h-1).
 * usage: hsdfs h [outfile]       | hsdfs h sample N seed | hsdfs h prefix N seed P
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

static const uint8_t MUL[4][4] = {{0,0,0,0},{0,1,2,3},{0,2,3,1},{0,3,1,2}};
static const uint8_t INV[4] = {0,1,3,2};

static int H, D;
static uint8_t V[65536][8];
static uint16_t SCAL[4][65536];
static uint8_t SMIN[65536];

static int vidx(const uint8_t *v) { int x = 0; for (int i = 0; i < 8; i++) x |= v[i] << (2*i); return x; }

static inline int sp(const uint8_t *a, const uint8_t *b) {
    int s = 0;
    for (int c = 0; c < D; c += 2) s ^= MUL[a[c]][b[c+1]] ^ MUL[a[c+1]][b[c]];
    return s;
}

static int nonsing(int k, uint8_t A[8][8]) {
    uint8_t B[8][8];
    memcpy(B, A, sizeof(B));
    for (int c = 0; c < k; c++) {
        int p = -1;
        for (int i = c; i < k; i++) if (B[i][c]) { p = i; break; }
        if (p < 0) return 0;
        if (p != c) for (int j = 0; j < k; j++) { uint8_t t = B[p][j]; B[p][j] = B[c][j]; B[c][j] = t; }
        uint8_t iv = INV[B[c][c]];
        for (int j = c; j < k; j++) B[c][j] = MUL[iv][B[c][j]];
        for (int i = c + 1; i < k; i++) if (B[i][c]) {
            uint8_t f = B[i][c];
            for (int j = c; j < k; j++) B[i][j] ^= MUL[f][B[c][j]];
        }
    }
    return 1;
}

static uint8_t M[8][8];
static int popc(int x) { return __builtin_popcount(x); }

static int minor_ok(int Rm, int Cm) {
    int rows[8], cols[8], nr = 0, nc = 0;
    for (int a = 0; a < H; a++) if (Rm >> a & 1) { rows[nr++] = 2*a; rows[nr++] = 2*a+1; }
    for (int b = 0; b < H; b++) if (Cm >> b & 1) { cols[nc++] = 2*b; cols[nc++] = 2*b+1; }
    uint8_t A[8][8];
    for (int i = 0; i < nr; i++) for (int j = 0; j < nc; j++) A[i][j] = M[rows[i]][cols[j]];
    return nonsing(nr, A);
}

static int party_ok(int p) {
    for (int Rm = 1 << p; Rm < (2 << p); Rm++) {
        int k = popc(Rm);
        for (int Cm = 1; Cm < (1 << H); Cm++) if (popc(Cm) == k)
            if (!minor_ok(Rm, Cm)) return 0;
    }
    return 1;
}

static int full_check_cnt(int *bad) {
    int nP = 2 * H;
    for (int S = 0; S < (1 << nP); S++) {
        if (popc(S) != H) continue;
        uint8_t A[8][8];
        int col = 0;
        for (int j = 0; j < nP; j++) if (S >> j & 1) {
            for (int t = 0; t < 2; t++) {
                int c = 2*j + t;
                for (int r = 0; r < D; r++) A[r][col] = (c < D) ? (r == c) : M[r][c - D];
                col++;
            }
        }
        if (!nonsing(D, A)) { if (!bad) return 0; (*bad)++; }
    }
    return bad ? (*bad == 0) : 1;
}
static int full_check(void) { return full_check_cnt(NULL); }

static int symp_ok(void) {
    for (int i = 0; i < D; i++) for (int j = i + 1; j < D; j++) {
        int want = (j == i + 1 && (i % 2 == 0));
        if (sp(M[i], M[j]) != want) return 0;
    }
    return 1;
}

static long nodes[16], nsol = 0, nstart = 0;
static FILE *fout = NULL;

static void dfs_party(int p, int *pool, int np) {
    if (p == H) {
        nodes[H]++;
        if (!symp_ok()) { fprintf(stderr, "BUG: not symplectic\n"); exit(2); }
        if (!full_check()) { fprintf(stderr, "BUG: full_check fails but party minors passed\n"); exit(2); }
        nsol++;
        if (fout) {
            for (int i = 0; i < D; i++) { for (int j = 0; j < D; j++) fprintf(fout, "%d", M[i][j]); fprintf(fout, " "); }
            fprintf(fout, "\n");
        }
        return;
    }
    int *np2 = malloc(sizeof(int) * (np + 1));
    for (int a = 0; a < np; a++) {
        int y0 = pool[a];
        if (!SMIN[y0]) continue;
        for (int b = 0; b < np; b++) {
            int y1 = pool[b];
            if (sp(V[y0], V[y1]) != 1) continue;
            /* canonical frame: y0 = min nonzero of plane; y1 = min of y1 + c*y0 */
            int ok = 1;
            for (int c = 1; c < 4 && ok; c++) if ((y1 ^ SCAL[c][y0]) < y1) ok = 0;
            if (!ok) continue;
            for (int s = 1; s < 4 && ok; s++) for (int t = 0; t < 4 && ok; t++) {
                int w = SCAL[s][y1] ^ SCAL[t][y0];
                if (w < y0) ok = 0;
            }
            for (int t = 1; t < 4 && ok; t++) if (SCAL[t][y0] < y0) ok = 0;
            if (!ok) continue;
            memcpy(M[2*p], V[y0], D);
            memcpy(M[2*p+1], V[y1], D);
            nodes[p]++;
            if (!party_ok(p)) continue;
            nodes[p+8]++;
            int n2 = 0;
            if (p + 1 < H) {
                for (int c = 0; c < np; c++) {
                    int z = pool[c];
                    if (sp(V[z], V[y0]) == 0 && sp(V[z], V[y1]) == 0) np2[n2++] = z;
                }
            }
            dfs_party(p + 1, np2, n2);
        }
    }
    free(np2);
}

static int tup[8];

static void start_party0(void) {
    for (int c = 0; c < H; c++) {
        M[0][2*c] = 1; M[0][2*c+1] = 0; M[1][2*c] = 0; M[1][2*c+1] = tup[c];
    }
    nstart++;
    int *pool = malloc(sizeof(int) * 65536);
    int np = 0;
    int lim = 1 << (2 * D);
    for (int v = 1; v < lim; v++)
        if (sp(V[v], M[0]) == 0 && sp(V[v], M[1]) == 0) pool[np++] = v;
    dfs_party(1, pool, np);
    free(pool);
}

static void rec_tuples(int c, int sum) {
    if (c == H) { if (sum == 1) start_party0(); return; }
    for (int d = 1; d < 4; d++) { tup[c] = d; rec_tuples(c + 1, sum ^ d); }
}

static uint64_t rs;
/* splitmix64: nonlinear over F_2 (plain xorshift is F_2-linear and can bias algebraic events) */
static inline uint64_t rnd(void) { uint64_t z = (rs += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return (z ^ (z >> 31)) >> 20; }

/* sample the first nr rows (nr even) of a uniform symplectic frame */
static int sample_rows(int nr) {
    for (int p = 0; p < nr / 2; p++) {
        int tries = 0, v;
        for (;;) {
            if (++tries > 2000000) return 0;
            v = (int)(rnd() & ((1 << (2*D)) - 1));
            if (!v) continue;
            int ok = 1;
            for (int q = 0; q < 2*p && ok; q++) if (sp(V[v], M[q]) != 0) ok = 0;
            if (ok) break;
        }
        memcpy(M[2*p], V[v], D);
        for (;;) {
            if (++tries > 2000000) return 0;
            int w = (int)(rnd() & ((1 << (2*D)) - 1));
            int ok = (sp(V[w], V[v]) == 1);
            for (int q = 0; q < 2*p && ok; q++) if (sp(V[w], M[q]) != 0) ok = 0;
            if (ok) { memcpy(M[2*p+1], V[w], D); break; }
        }
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
    H = atoi(argv[1]); D = 2 * H;
    for (int x = 0; x < 65536; x++) for (int i = 0; i < 8; i++) V[x][i] = (x >> (2*i)) & 3;
    for (int s = 0; s < 4; s++) for (int x = 0; x < 65536; x++) {
        uint8_t w[8]; for (int i = 0; i < 8; i++) w[i] = MUL[s][V[x][i]];
        SCAL[s][x] = vidx(w);
    }
    for (int x = 0; x < (1 << (2*D)); x++) SMIN[x] = (x != 0 && x <= SCAL[2][x] && x <= SCAL[3][x]);
    if (argc >= 5 && (!strcmp(argv[2], "sample") || !strcmp(argv[2], "prefix"))) {
        int prefix = !strcmp(argv[2], "prefix");
        long N = atol(argv[3]); rs = 88172645463325252ULL ^ (uint64_t)atol(argv[4]) * 0x9E3779B97F4A7C15ULL;
        int P = prefix ? atoi(argv[5]) : H - 1;
        for (int i = 0; i < 20; i++) rnd();
        long good = 0, got = 0;
        for (long s = 0; s < N; s++) {
            if (!sample_rows(2 * (P + 1))) continue;
            got++;
            int pass = 1;
            if (prefix) { for (int q = 0; q <= P && pass; q++) if (!party_ok(q)) pass = 0; }
            else pass = full_check();
            if (pass) good++;
        }
        long double cnt = 1;
        for (int p = 0; p <= P; p++) { long double a = 1, b = 1; for (int i = 0; i < D - 2*p; i++) a *= 4; for (int i = 0; i < D - 2*p - 1; i++) b *= 4; cnt *= (a - 1) * b; }
        long double sc = 1; for (int i = 0; i < H + P; i++) sc *= 60;
        long double pr = (long double)good / got;
        printf("%s h=%d P=%d N=%ld good=%ld frac=%.6Lf est normalised nodes = %.2Lf (+- %.2Lf)\n", argv[2], H, P, got, good, pr, pr * cnt / sc,
               good ? sqrtl((long double)good) / got * cnt / sc : 0.0L);
        return 0;
    }
    if (argc >= 3 && strcmp(argv[2], "-")) fout = fopen(argv[2], "w");
    rec_tuples(0, 0);
    printf("S-type h=%d starts=%ld N_norm=%ld (all solutions = N_norm*60^%d)\n", H, nstart, nsol, 2*H-1);
    for (int p = 1; p < H; p++) printf("  party %d: candidates %ld passed-minors %ld\n", p, nodes[p], nodes[p+8]);
    if (fout) fclose(fout);
    return 0;
}
