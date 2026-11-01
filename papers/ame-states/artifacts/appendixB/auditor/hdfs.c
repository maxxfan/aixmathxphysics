/* hdfs.c -- auditor's independent exhaustive search for H-type AME(n,4), n = 2h parties.
 *
 * Model: stabilizer = rowspace of G = [I_{2h} | M] over F_4, M is (2h x 2h), M M^dagger = I  (Hermitian self-dual),
 * party j = coordinates 2j,2j+1 of the 4h-long codeword (columns of G: 0..2h-1 are the identity part).
 * AME  <=>  every party block-minor M[R,C], |R|=|C|=k, is nonsingular  (equivalently all C(2h,h) column sets of
 * G have rank 2h; the final full check does the latter directly).
 *
 * Normal form (all steps are equivalences of the problem):
 *   - right-multiplying M by a block-diagonal unitary (local U(2,4), order 18 per column party) fixes the blocks of
 *     row party 0 up to H_c = X_c X_c^dagger (10 invertible 2x2 Hermitian matrices); take a fixed representative X_H.
 *     Condition: sum_c H_c = I (row party 0 orthonormal).
 *   - left-multiplying by block-diag unitary on row parties 1..h-1: each party's frame (r0,r1) is canonical:
 *     r0 minimal among its F_4^* multiples, r1 minimal among its multiples, r0 < r1  (the 6 norm-1 vectors of a
 *     nondegenerate Hermitian plane are exactly the multiples of r0 and of r1).
 *   Hence #all solutions = N_norm * 18^h * 18^(h-1).   Check: h=2 -> N_norm = 34992/5832 = 6.
 * usage: hdfs h [outfile]            exhaustive DFS
 *        hdfs h sample N seed        random unitary M, estimate #solutions (cross-check of the normalisation)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

static const uint8_t MUL[4][4] = {{0,0,0,0},{0,1,2,3},{0,2,3,1},{0,3,1,2}};
static const uint8_t INV[4] = {0,1,3,2};
static const uint8_t CONJ[4] = {0,1,3,2};

static int H;           /* h */
static int D;           /* 2h = dimension of M, also number of parties */
static uint8_t V[65536][8];
static uint8_t SM[4][65536];  /* not used: scalar multiples via idx tables */
static uint16_t SCAL[4][65536];
static uint8_t NORM1[65536];
static uint8_t SMIN[65536];

static int vidx(const uint8_t *v) { int x = 0; for (int i = 0; i < 8; i++) x |= v[i] << (2*i); return x; }

static inline int ip(const uint8_t *a, const uint8_t *b) {
    int s = 0;
    for (int i = 0; i < D; i++) s ^= MUL[a[i]][CONJ[b[i]]];
    return s;
}

/* rank of k x k matrix (k<=8) */
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

/* party-block minor M[Rmask, Cmask] nonsingular? masks over h parties */
static int minor_ok(int Rm, int Cm) {
    int rows[8], cols[8], nr = 0, nc = 0;
    for (int a = 0; a < H; a++) if (Rm >> a & 1) { rows[nr++] = 2*a; rows[nr++] = 2*a+1; }
    for (int b = 0; b < H; b++) if (Cm >> b & 1) { cols[nc++] = 2*b; cols[nc++] = 2*b+1; }
    uint8_t A[8][8];
    for (int i = 0; i < nr; i++) for (int j = 0; j < nc; j++) A[i][j] = M[rows[i]][cols[j]];
    return nonsing(nr, A);
}

static int popc(int x) { return __builtin_popcount(x); }

/* check all minors with row set containing party p, row set within 0..p */
static int KMAX = 99, STOPAT = 99;
static int party_ok(int p) {
    for (int Rm = 1 << p; Rm < (2 << p); Rm++) {
        if (!(Rm >> p & 1)) continue;
        int k = popc(Rm);
        if (k > KMAX) continue;
        for (int Cm = 1; Cm < (1 << H); Cm++) if (popc(Cm) == k)
            if (!minor_ok(Rm, Cm)) return 0;
    }
    return 1;
}

/* full direct check: every set S of h parties among 2h: the 2h columns of G=[I|M] have rank 2h.
   If bad != NULL, count failing sets S instead of stopping at the first. */
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
static long badhist[80];

static int unitary_ok(void) {
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
        int s = 0;
        for (int c = 0; c < D; c++) s ^= MUL[M[i][c]][CONJ[M[j][c]]];
        if (s != (i == j)) return 0;
    }
    return 1;
}

/* ---------------- DFS ---------------- */
static long nodes[16], nsol = 0;
static FILE *fout = NULL;

typedef struct { int n; int *idx; } Pool;

static int *pool_buf[16];

static void dfs_party(int p, int *pool, int np) {
    if (p == H) {
        nodes[H]++;
        if (!unitary_ok()) { fprintf(stderr, "BUG: not unitary\n"); exit(2); }
        if (!full_check()) { fprintf(stderr, "BUG: full_check fails but party minors passed\n"); exit(2); }
        nsol++;
        if (fout) {
            for (int i = 0; i < D; i++) { for (int j = 0; j < D; j++) fprintf(fout, "%d", M[i][j]); fprintf(fout, " "); }
            fprintf(fout, "\n");
        }
        return;
    }
    int *tmp = malloc(sizeof(int) * (np + 1));
    int *np2 = malloc(sizeof(int) * (np + 1));
    for (int a = 0; a < np; a++) {
        int y0 = pool[a];
        /* candidates for y1: in pool, orthogonal to y0, index > y0 */
        int nt = 0;
        for (int b = a + 1; b < np; b++) {
            int y1 = pool[b];
            if (ip(V[y0], V[y1]) == 0) tmp[nt++] = y1;
        }
        for (int t = 0; t < nt; t++) {
            int y1 = tmp[t];
            memcpy(M[2*p], V[y0], D);
            memcpy(M[2*p+1], V[y1], D);
            nodes[p]++;
            if (!party_ok(p)) continue;
            nodes[p+8]++;
            if (p == STOPAT) {
                if (p == H - 2) {
                    /* complete the last party (unique canonical frame) and histogram the number of failing party sets */
                    int zl[4096], nz = 0;
                    for (int b = 0; b < np; b++) {
                        int v = pool[b];
                        if (ip(V[v], M[2*p]) == 0 && ip(V[v], M[2*p+1]) == 0 && nz < 4096) zl[nz++] = v;
                    }
                    int found = 0;
                    for (int a2 = 0; a2 < nz && !found; a2++) for (int b2 = a2 + 1; b2 < nz && !found; b2++)
                        if (ip(V[zl[a2]], V[zl[b2]]) == 0) {
                            memcpy(M[2*(p+1)], V[zl[a2]], D); memcpy(M[2*(p+1)+1], V[zl[b2]], D);
                            int bad = 0; full_check_cnt(&bad);
                            badhist[bad]++;
                            found = 1;
                        }
                    if (!found) badhist[79]++;
                }
                continue;
            }
            /* new pool: pool members orthogonal to y0 and y1 */
            int n2 = 0;
            if (p + 1 < H) {
                for (int b = 0; b < np; b++) {
                    int z = pool[b];
                    if (ip(V[z], V[y0]) == 0 && ip(V[z], V[y1]) == 0) np2[n2++] = z;
                }
            }
            dfs_party(p + 1, np2, n2);
        }
    }
    free(tmp); free(np2);
}

/* 2x2 Hermitian invertible reps */
static int nH = 0;
static uint8_t HX[10][2][2];   /* rep X with X X^dagger = H */
static uint8_t HH[10][2][2];

static void build_H(void) {
    /* enumerate all 2x2 matrices X over F_4; H = X X^dagger; keep invertible ones, first rep per H */
    int seen[256]; memset(seen, -1, sizeof(seen));
    for (int x = 0; x < 256; x++) {
        uint8_t X[2][2] = {{x & 3, (x >> 2) & 3}, {(x >> 4) & 3, (x >> 6) & 3}};
        uint8_t Hm[2][2];
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
            int s = 0;
            for (int c = 0; c < 2; c++) s ^= MUL[X[i][c]][CONJ[X[j][c]]];
            Hm[i][j] = s;
        }
        int det = MUL[Hm[0][0]][Hm[1][1]] ^ MUL[Hm[0][1]][Hm[1][0]];
        if (!det) continue;
        int key = Hm[0][0] | (Hm[0][1] << 2) | (Hm[1][0] << 4) | (Hm[1][1] << 6);
        if (seen[key] >= 0) continue;
        seen[key] = nH;
        memcpy(HX[nH], X, 4); memcpy(HH[nH], Hm, 4);
        nH++;
    }
}

static int tup[8];
static long nstart = 0;

static void start_party0(void) {
    for (int c = 0; c < H; c++) {
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) M[i][2*c + j] = HX[tup[c]][i][j];
    }
    nstart++;
    /* pool: norm-1, scalar-min, orthogonal to rows 0,1 (rows >= 2 zero for now: ip uses only first D entries) */
    int *pool = malloc(sizeof(int) * 65536);
    int np = 0;
    int lim = 1 << (2 * D);
    for (int v = 0; v < lim; v++) {
        if (!SMIN[v] || !NORM1[v]) continue;
        if (ip(V[v], M[0]) == 0 && ip(V[v], M[1]) == 0) pool[np++] = v;
    }
    dfs_party(1, pool, np);
    free(pool);
}

static void rec_tuples(int c, int sum[2][2]) {
    if (c == H) {
        /* sum must be I */
        if (sum[0][0] == 1 && sum[0][1] == 0 && sum[1][0] == 0 && sum[1][1] == 1) start_party0();
        return;
    }
    for (int h = 0; h < nH; h++) {
        int s2[2][2];
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) s2[i][j] = sum[i][j] ^ HH[h][i][j];
        tup[c] = h;
        if (c == H - 1) { if (!(s2[0][0] == 1 && s2[0][1] == 0 && s2[1][0] == 0 && s2[1][1] == 1)) continue; }
        rec_tuples(c + 1, s2);
    }
}

/* ---------------- sampling mode ---------------- */
static uint64_t rs;
/* splitmix64: nonlinear over F_2 (plain xorshift is F_2-linear and can bias algebraic events) */
static inline uint64_t rnd(void) { uint64_t z = (rs += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return (z ^ (z >> 31)) >> 20; }

static int sample_unitary(void) {
    /* rows sequentially: random norm-1 vector orthogonal to previous rows (rejection) */
    for (int r = 0; r < D; r++) {
        int tries = 0;
        for (;;) {
            if (++tries > 200000) return 0;
            int v = (int)(rnd() & ((1 << (2*D)) - 1));
            if (!NORM1[v]) continue;
            int ok = 1;
            for (int q = 0; q < r; q++) if (ip(V[v], M[q]) != 0) { ok = 0; break; }
            if (!ok) continue;
            memcpy(M[r], V[v], D);
            break;
        }
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: hdfs h [outfile] | hdfs h sample N seed\n"); return 1; }
    H = atoi(argv[1]); D = 2 * H;
    for (int x = 0; x < 65536; x++) for (int i = 0; i < 8; i++) V[x][i] = (x >> (2*i)) & 3;
    for (int s = 0; s < 4; s++) for (int x = 0; x < 65536; x++) {
        uint8_t w[8]; for (int i = 0; i < 8; i++) w[i] = MUL[s][V[x][i]];
        SCAL[s][x] = vidx(w);
    }
    for (int x = 0; x < (1 << (2*D)); x++) {
        NORM1[x] = (ip(V[x], V[x]) == 1);
        SMIN[x] = (x != 0 && x <= SCAL[2][x] && x <= SCAL[3][x]);
    }
    if (argc >= 6 && !strcmp(argv[2], "prefix")) {
        /* hdfs h prefix N seed P : sample orthonormal prefixes of rows 0..2P+1, test party_ok(0..P); estimate normalised node count */
        long N = atol(argv[3]); rs = 88172645463325252ULL ^ (uint64_t)atol(argv[4]) * 0x9E3779B97F4A7C15ULL;
        int P = atoi(argv[5]);
        if (argc >= 7) KMAX = atoi(argv[6]);
        for (int i = 0; i < 20; i++) rnd();
        long good = 0, got = 0;
        int nr = 2 * (P + 1);
        for (long s = 0; s < N; s++) {
            int okk = 1;
            for (int r = 0; r < nr && okk; r++) {
                int tries = 0;
                for (;;) {
                    if (++tries > 200000) { okk = 0; break; }
                    int v = (int)(rnd() & ((1 << (2*D)) - 1));
                    if (!NORM1[v]) continue;
                    int ok = 1;
                    for (int q = 0; q < r; q++) if (ip(V[v], M[q]) != 0) { ok = 0; break; }
                    if (!ok) continue;
                    memcpy(M[r], V[v], D);
                    break;
                }
            }
            if (!okk) continue;
            got++;
            int pass = 1;
            for (int q = 0; q <= P && pass; q++) if (!party_ok(q)) pass = 0;
            if (pass) good++;
        }
        long double cnt = 1;
        for (int r = 0; r < nr; r++) { int m = D - r; long double t = 1; for (int i = 0; i < m - 1; i++) t *= 2; t *= ((1L << m) - ((m % 2) ? -1 : 1)); cnt *= t; }
        long double sc = 1; for (int i = 0; i < H + P; i++) sc *= 18;
        long double p = (long double)good / got;
        printf("prefix h=%d P=%d N=%ld good=%ld frac=%.6Lf  est normalised nodes = %.1Lf (+- %.1Lf)\n", H, P, got, good, p, p * cnt / sc,
               (good ? sqrtl((long double)good) / got * cnt / sc : 0));
        return 0;
    }
    if (argc >= 5 && !strcmp(argv[2], "sample")) {
        long N = atol(argv[3]); rs = 88172645463325252ULL ^ (uint64_t)atol(argv[4]) * 0x9E3779B97F4A7C15ULL;
        for (int i = 0; i < 20; i++) rnd();
        long good = 0, got = 0;
        for (long s = 0; s < N; s++) {
            if (!sample_unitary()) continue;
            got++;
            if (full_check()) good++;
        }
        /* |GU(D,2)| = 2^(D(D-1)/2) prod_{i=1}^D (2^i - (-1)^i) */
        long double gu = 1; for (int i = 0; i < D*(D-1)/2; i++) gu *= 2;
        for (int i = 1; i <= D; i++) gu *= ((1L << i) - ((i % 2) ? -1 : 1));
        long double sc = 1; for (int i = 0; i < H + (H - 1); i++) sc *= 18;
        long double p = (long double)good / got;
        printf("sample h=%d N=%ld good=%ld frac=%.6Lf  |GU|=%.6Le  est N_norm = frac*|GU|/18^(2h-1) = %.2Lf (+- %.2Lf)\n",
               H, got, good, p, gu, p * gu / sc, (good ? sqrtl((long double)good) / got * gu / sc : 0));
        return 0;
    }
    if (argc >= 3 && strcmp(argv[2], "-")) fout = fopen(argv[2], "w");
    if (argc >= 5) { KMAX = atoi(argv[3]); STOPAT = atoi(argv[4]); }
    build_H();
    fprintf(stderr, "invertible Hermitian 2x2: %d\n", nH);
    int z[2][2] = {{0,0},{0,0}};
    rec_tuples(0, z);
    printf("h=%d starts=%ld N_norm=%ld  (all solutions = N_norm*18^%d)\n", H, nstart, nsol, 2*H-1);
    for (int p = 1; p < H; p++) printf("  party %d: candidates %ld passed-minors %ld\n", p, nodes[p], nodes[p+8]);
    for (int b = 0; b < 80; b++) if (badhist[b]) printf("  bad-set histogram: %d failing party sets (of 70) : %ld nodes\n", b, badhist[b]);
    if (fout) fclose(fout);
    return 0;
}
