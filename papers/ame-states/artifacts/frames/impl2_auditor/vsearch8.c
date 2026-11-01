/* Auditor phase 2: independent n=8 search in the V-picture (classify #315), own code.
   Given V_0 (4 basis vectors, 16-bit, party p = nibble p of the 16-bit word, W_4..W_7 = nibbles 0..3):
   count ordered (V_1,V_2,V_3): mutually orthogonal nondegenerate 4-dim subspaces of B (16-dim) with B = sum V_j and
     every nonzero element of every V_j has all 4 nibbles nonzero          (V_j cap W_L = 0, |L| = 3)
     every nonzero element of V_i + V_j has at least 3 nonzero nibbles     (|L| = 2), all pairs i<j among 0..3
     (|J| = 3 conditions follow by symplectic duality; V_3 conditions are checked directly, including its pairs).
   usage: vsearch8 b0 b1 b2 b3 [mode]   mode 0 = full (V_3 check), 1 = stop at V_2 count (relaxed)
   prints: n_V1, n_V2pairs, n_solutions ; writes solutions to stdout lines "SOL ..."
   cc -O2 -o vsearch8 vsearch8.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline int sw16(int v) { return ((v >> 2) & 0x3333) | ((v & 0x3333) << 2); }
static inline int om(int u, int v) { return __builtin_popcount(u & sw16(v)) & 1; }
static inline int wt(int x) { return __builtin_popcount((x | (x >> 1) | (x >> 2) | (x >> 3)) & 0x1111); }

static int nonsing4(const int *rows4) { /* rows: 4-bit */
    int r[4] = {rows4[0], rows4[1], rows4[2], rows4[3]}, rank = 0;
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
static int gram_nondeg(const int *v) {
    int g[4];
    for (int i = 0; i < 4; i++) { g[i] = 0; for (int j = 0; j < 4; j++) if (om(v[i], v[j])) g[i] |= 1 << j; }
    return nonsing4(g);
}

static int V0[16];      /* all 16 elements of V_0 */
static int V1[16], V2[16];
static unsigned char in_c[65536];

/* enumerate 4-dim subspaces entirely inside the candidate set (in_c flags, 0 excluded) via lexicographically least bases.
   cand: sorted ascending list of candidate vectors. callback receives basis[4]. */
typedef void (*cb_t)(const int *basis, void *ctx);

static void enum_sub(const int *cand, int nc, cb_t cb, void *ctx) {
    int basis[4];
    int span[16];
    for (int i1 = 0; i1 < nc; i1++) {
        int x1 = cand[i1];
        basis[0] = x1; span[0] = 0; span[1] = x1;
        for (int i2 = i1 + 1; i2 < nc; i2++) {
            int x2 = cand[i2];
            /* x2 must be min of coset x2+{0,x1}, and x2^x1 in candidates */
            if ((x2 ^ x1) < x2 || !in_c[x2 ^ x1]) continue;
            span[2] = x2; span[3] = x2 ^ x1;
            basis[1] = x2;
            for (int i3 = i2 + 1; i3 < nc; i3++) {
                int x3 = cand[i3];
                int ok = 1;
                for (int s = 0; s < 4; s++) { int y = x3 ^ span[s]; if (y < x3 || !in_c[y]) { ok = 0; break; } }
                if (!ok) continue;
                for (int s = 0; s < 4; s++) span[4 + s] = x3 ^ span[s];
                basis[2] = x3;
                for (int i4 = i3 + 1; i4 < nc; i4++) {
                    int x4 = cand[i4];
                    int ok4 = 1;
                    for (int s = 0; s < 8; s++) { int y = x4 ^ span[s]; if (y < x4 || !in_c[y]) { ok4 = 0; break; } }
                    if (!ok4) continue;
                    basis[3] = x4;
                    cb(basis, ctx);
                }
            }
        }
    }
}

static long long n_v1 = 0, n_v2 = 0, n_sol = 0, n_v1_with_v2 = 0;
static int mode = 0;
static long long n_v3_single = 0, n_v3_pair = 0, n_dual_viol = 0, n_pair_only = 0;
static int tS = 4, tP = 3;          /* mode 0: single-element weight 4, pair weight >= 3; mode 2 (relaxed control): 3 and 2 */
static long long maxv1 = 1LL << 60;  /* stop after this many V_1 (relaxed control only) */
static int Rlist[4096], nR;
static int B0[4];

typedef struct { int dummy; } ctx_t;

static int cur_b1[4];
static int R2[256], nR2;
static int found_v2_for_this_v1;

static void expand(const int *b, int *out) { /* span of 4 basis vectors -> 16 elements */
    for (int m = 0; m < 16; m++) { int x = 0; for (int i = 0; i < 4; i++) if ((m >> i) & 1) x ^= b[i]; out[m] = x; }
}

static void cb_v2(const int *b2, void *ctx) {
    (void)ctx;
    if (!gram_nondeg(b2)) return;
    n_v2++;
    found_v2_for_this_v1 = 1;
    if (mode == 1) return;
    expand(b2, V2);
    /* V_3 = R2 cap V_2^perp */
    int v3[16], n3 = 0;
    for (int i = 0; i < nR2; i++) {
        int x = R2[i], ok = 1;
        for (int k = 0; k < 4; k++) if (om(x, b2[k])) { ok = 0; break; }
        if (ok) v3[n3++] = x;
    }
    if (n3 != 16) { fprintf(stderr, "V3 size %d\n", n3); exit(2); }
    {   /* diagnostics: does V_3 pass its single-element weights, and then its pair conditions (duality says: automatically)? */
        int s_ok = 1, p_ok = 1;
        for (int i = 0; i < 16; i++) {
            if (v3[i] == 0) continue;
            if (wt(v3[i]) < tS) s_ok = 0;
            for (int j = 0; j < 16; j++)
                if (wt(v3[i] ^ V0[j]) < tP || wt(v3[i] ^ V1[j]) < tP || wt(v3[i] ^ V2[j]) < tP) p_ok = 0;
        }
        if (s_ok) n_v3_single++;
        if (p_ok) n_v3_pair++;
        if (s_ok && !p_ok) n_dual_viol++;
        if (p_ok && !s_ok) {
            n_pair_only++;
            if (n_pair_only <= 30)
                printf("SOL V0=%d,%d,%d,%d V1=%d,%d,%d,%d V2=%d,%d,%d,%d\n", B0[0], B0[1], B0[2], B0[3],
                       cur_b1[0], cur_b1[1], cur_b1[2], cur_b1[3], b2[0], b2[1], b2[2], b2[3]);
        }
    }
    for (int i = 0; i < 16; i++) {
        if (v3[i] == 0) continue;
        if (wt(v3[i]) < tS) return;
        for (int j = 0; j < 16; j++) {
            if (wt(v3[i] ^ V0[j]) < tP || wt(v3[i] ^ V1[j]) < tP || wt(v3[i] ^ V2[j]) < tP) return;
        }
    }
    n_sol++;
    printf("SOL V0=%d,%d,%d,%d V1=%d,%d,%d,%d V2=%d,%d,%d,%d\n", B0[0], B0[1], B0[2], B0[3],
           cur_b1[0], cur_b1[1], cur_b1[2], cur_b1[3], b2[0], b2[1], b2[2], b2[3]);
    fflush(stdout);
}

static void cb_v1(const int *b1, void *ctx) {
    (void)ctx;
    if (!gram_nondeg(b1)) return;
    n_v1++;
    if (n_v1 > maxv1) {
        printf("RESULT (capped at %lld V_1) n_V1=%lld n_V1_with_V2=%lld n_V2=%lld n_SOL=%lld\n", maxv1, n_v1 - 1, n_v1_with_v2, n_v2, n_sol);
        exit(0);
    }
    memcpy(cur_b1, b1, sizeof cur_b1);
    expand(b1, V1);
    /* R2 = R cap V_1^perp */
    nR2 = 0;
    for (int i = 0; i < nR; i++) {
        int x = Rlist[i], ok = 1;
        for (int k = 0; k < 4; k++) if (om(x, b1[k])) { ok = 0; break; }
        if (ok) R2[nR2++] = x;
    }
    if (nR2 != 256) { fprintf(stderr, "R2 size %d\n", nR2); exit(2); }
    /* candidates for V_2 */
    static unsigned char saved[65536];
    int cand2[256], nc2 = 0;
    /* reuse in_c: save current in_c flags of R2 elements */
    memset(saved, 0, sizeof saved);
    for (int i = 0; i < nR2; i++) saved[R2[i]] = 1;
    static unsigned char c1flag[65536];
    memcpy(c1flag, in_c, sizeof c1flag);
    memset(in_c, 0, sizeof in_c);
    for (int i = 0; i < nR2; i++) {
        int x = R2[i];
        if (x == 0 || wt(x) < tS) continue;
        int ok = 1;
        for (int j = 0; j < 16 && ok; j++) if (wt(x ^ V0[j]) < tP || wt(x ^ V1[j]) < tP) ok = 0;
        if (ok) { in_c[x] = 1; cand2[nc2++] = x; }
    }
    /* sort ascending */
    for (int i = 1; i < nc2; i++) { int v = cand2[i], j = i - 1; while (j >= 0 && cand2[j] > v) { cand2[j + 1] = cand2[j]; j--; } cand2[j + 1] = v; }
    found_v2_for_this_v1 = 0;
    enum_sub(cand2, nc2, cb_v2, NULL);
    if (found_v2_for_this_v1) n_v1_with_v2++;
    memcpy(in_c, c1flag, sizeof in_c);
}

int main(int argc, char **argv) {
    if (argc < 5) return 1;
    for (int i = 0; i < 4; i++) B0[i] = atoi(argv[1 + i]);
    if (argc > 5) mode = atoi(argv[5]);
    if (mode == 2) { tS = 3; tP = 2; }
    if (argc > 6) maxv1 = atoll(argv[6]);
    expand(B0, V0);
    /* sanity: V_0 isotropy of the setup is not required; check nondegenerate and weights */
    if (!gram_nondeg(B0)) { fprintf(stderr, "V0 degenerate\n"); return 2; }
    for (int i = 1; i < 16; i++) if (wt(V0[i]) != 4) { fprintf(stderr, "V0 weight %d\n", wt(V0[i])); return 2; }
    nR = 0;
    for (int x = 0; x < 65536; x++) {
        int ok = 1;
        for (int k = 0; k < 4; k++) if (om(x, B0[k])) { ok = 0; break; }
        if (ok) Rlist[nR++] = x;
    }
    if (nR != 4096) { fprintf(stderr, "R size %d\n", nR); return 2; }
    int cand1[4096], nc1 = 0;
    memset(in_c, 0, sizeof in_c);
    for (int i = 0; i < nR; i++) {
        int x = Rlist[i];
        if (x == 0 || wt(x) < tS) continue;
        int ok = 1;
        for (int j = 0; j < 16 && ok; j++) if (wt(x ^ V0[j]) < tP) ok = 0;
        if (ok) { in_c[x] = 1; cand1[nc1++] = x; }
    }
    /* Rlist is generated in increasing x already; cand1 ascending */
    fprintf(stderr, "cand1 = %d\n", nc1);
    enum_sub(cand1, nc1, cb_v1, NULL);
    printf("RESULT n_V1=%lld n_V1_with_V2=%lld n_V2=%lld n_SOL=%lld | V3 singles ok %lld, V3 pairs ok %lld, singles-ok-but-pairs-fail %lld, pairs-ok-but-singles-fail %lld\n",
           n_v1, n_v1_with_v2, n_v2, n_sol, n_v3_single, n_v3_pair, n_dual_viol, n_pair_only);
    return 0;
}
