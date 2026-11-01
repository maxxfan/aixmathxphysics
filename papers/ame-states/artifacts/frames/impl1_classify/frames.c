/* frames.c : exhaustive search for frames of an F_2-linear stabilizer AME(2m,4) state (m = 2,3,4).
 *
 * Space B = W_0 + ... + W_{m-1}, W_l = F_2^4 with form J (coordinates x1,x2,z1,z2).  A frame is an ordered list of
 * mutually orthogonal nondegenerate 4-dim subspaces V_0..V_{m-1} of B, each a graph over W_0:
 *     V_i = { (w, th_{i,1} w, ..., th_{i,m-1} w) : w in F_2^4 },  th_{i,l} in GL(4,2)   (th_{i,0} = I)
 * with   V_J ∩ W_L = 0  for all J,L (subsets of {0..m-1}) with |J|+|L| = t-1      (t = min party-weight; AME: t = m+1).
 * Level 0 is read from a list of gauge-class representatives (v0reps.py); levels 1..m-1 are enumerated completely:
 * the orthogonality equations  sum_{l=0}^{m-1} th_{i,l}^T J th_{j,l} = 0  (i<j, th_{.,0} = I) are linear in th_{j,.};
 * their affine solution space (dim f) is enumerated in Gray-code order, every point is tested.
 *
 * usage: frames m t v0file classindex [maxsol] [outfile]
 *   classindex = -1 : all classes in the file.  Prints node counts per level (nodes[j] = number of valid V_0..V_j).
 * Packed matrices: entry (r,c) = bit 4r+c of a uint16.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef uint16_t u16;
typedef uint64_t u64;

static int M_, T_;
static u16 JM;
static u16 FORM[65536];   /* th^T J th */
static unsigned char GL[65536];

static u16 mul4(u16 a, u16 b) {
    u16 c = 0;
    for (int r = 0; r < 4; r++) {
        int row = (a >> (4 * r)) & 15, acc = 0;
        for (int k = 0; k < 4; k++) if ((row >> k) & 1) acc ^= (b >> (4 * k)) & 15;
        c |= acc << (4 * r);
    }
    return c;
}
static u16 tr4(u16 a) {
    u16 c = 0;
    for (int r = 0; r < 4; r++) for (int k = 0; k < 4; k++) if ((a >> (4 * r + k)) & 1) c |= 1 << (4 * k + r);
    return c;
}
static int rank4(u16 a) {
    int rows[4], rk = 0;
    for (int r = 0; r < 4; r++) rows[r] = (a >> (4 * r)) & 15;
    for (int bit = 0; bit < 4; bit++) {
        int piv = -1;
        for (int i = rk; i < 4; i++) if ((rows[i] >> bit) & 1) { piv = i; break; }
        if (piv < 0) continue;
        int t = rows[rk]; rows[rk] = rows[piv]; rows[piv] = t;
        for (int i = 0; i < 4; i++) if (i != rk && ((rows[i] >> bit) & 1)) rows[i] ^= rows[rk];
        rk++;
    }
    return rk;
}
static void init_tables(void) {
    JM = 0;
    int rc[4][2] = {{0, 2}, {1, 3}, {2, 0}, {3, 1}};
    for (int i = 0; i < 4; i++) JM |= 1 << (4 * rc[i][0] + rc[i][1]);
    for (int a = 0; a < 65536; a++) {
        FORM[a] = mul4(mul4(tr4((u16)a), JM), (u16)a);
        GL[a] = rank4((u16)a) == 4;
    }
}

/* frame storage: th[i][l] for member i, party l (th[i][0] = identity) */
static u16 th[8][8];
static u16 I4;
static u16 vec[8][4];         /* basis vectors of V_i in F_2^{4m}, party l in bits 4l..4l+3 */
static long long nodes[8];
static long long nsol, maxsol = -1;
static FILE *fout;

static void make_vecs(int i) {
    for (int w = 0; w < 4; w++) {
        u64 v = 0;
        for (int l = 0; l < M_; l++) {
            int col = 0;
            for (int s = 0; s < 4; s++) col |= ((th[i][l] >> (4 * s + w)) & 1) << s;
            v |= (u64)col << (4 * l);
        }
        vec[i][w] = (u16)v;
    }
}

/* rank of up to 16 vectors masked by mask */
static int rank_vecs(const u16 *v, int cnt, u16 mask) {
    u16 basis[16]; int nb = 0;
    for (int k = 0; k < cnt; k++) {
        u16 x = v[k] & mask;
        for (int b = 0; b < nb; b++) { u16 y = x ^ basis[b]; if (y < x) x = y; }
        if (x) {
            basis[nb++] = x;
            /* keep basis sorted descending so that the min-trick works: insert-sort */
            for (int b = nb - 1; b > 0 && basis[b] > basis[b - 1]; b--) { u16 t = basis[b]; basis[b] = basis[b - 1]; basis[b - 1] = t; }
        } else return -1;  /* dependent => injectivity fails */
    }
    return nb;
}

/* conditions for level j: all (Jmask ∋ j, Jmask ⊂ {0..j}), Lmask with popcount(J)+popcount(L) = T_-1, L nonempty */
typedef struct { int jm, lm; } cond_t;
static cond_t conds[8][512];
static int ncond[8];
static void build_conds(void) {
    for (int j = 0; j < M_; j++) {
        ncond[j] = 0;
        for (int jm = 1; jm < (1 << (j + 1)); jm++) {
            if (!((jm >> j) & 1)) continue;
            int pj = __builtin_popcount(jm);
            for (int lm = 1; lm < (1 << M_); lm++) {
                if (pj + __builtin_popcount(lm) != T_ - 1) continue;
                conds[j][ncond[j]].jm = jm; conds[j][ncond[j]].lm = lm; ncond[j]++;
            }
        }
    }
}
static int check_conds(int j) {
    make_vecs(j);
    for (int c = 0; c < ncond[j]; c++) {
        int jm = conds[j][c].jm, lm = conds[j][c].lm;
        u16 v[16]; int cnt = 0;
        for (int i = 0; i <= j; i++) if ((jm >> i) & 1) for (int w = 0; w < 4; w++) v[cnt++] = vec[i][w];
        u16 mask = 0;
        for (int l = 0; l < M_; l++) if (!((lm >> l) & 1)) mask |= 15 << (4 * l);
        if (rank_vecs(v, cnt, mask) < 0) return 0;   /* V_J -> W_{L^c} not injective */
    }
    return 1;
}

static void print_frame(void) {
    for (int i = 0; i < M_; i++) {
        for (int l = 1; l < M_; l++) fprintf(fout, "%u ", th[i][l]);
        fprintf(fout, "| ");
    }
    fprintf(fout, "\n");
}

static void level(int j);

/* solve the linear system for level j and enumerate */
static void level(int j) {
    if (j == M_) {
        nsol++;
        if (fout) print_frame();
        return;
    }
    int nu = 16 * (M_ - 1);          /* unknowns */
    int neq = 16 * j;
    u64 eq[64]; int ne = 0;
    const u64 RHS = 1ULL << 62;
    for (int i = 0; i < j; i++) {
        for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) {
            u64 row = 0;
            for (int l = 1; l < M_; l++) {
                u16 C = mul4(tr4(th[i][l]), JM);   /* (T_l^T J)[a][s] */
                for (int s = 0; s < 4; s++) if ((C >> (4 * a + s)) & 1) row |= 1ULL << (16 * (l - 1) + 4 * s + b);
            }
            if ((JM >> (4 * a + b)) & 1) row |= RHS;
            eq[ne++] = row;
        }
    }
    (void)neq;
    /* Gaussian elimination (reduced echelon) */
    int piv[64], npiv = 0;
    int isp[64]; memset(isp, 0, sizeof isp);
    int r = 0;
    for (int c = 0; c < nu && r < ne; c++) {
        int p = -1;
        for (int i = r; i < ne; i++) if ((eq[i] >> c) & 1) { p = i; break; }
        if (p < 0) continue;
        u64 t = eq[r]; eq[r] = eq[p]; eq[p] = t;
        for (int i = 0; i < ne; i++) if (i != r && ((eq[i] >> c) & 1)) eq[i] ^= eq[r];
        piv[npiv++] = c; isp[c] = 1; r++;
    }
    for (int i = r; i < ne; i++) if (eq[i] & RHS) return;   /* inconsistent */
    u64 xp = 0;
    for (int k = 0; k < npiv; k++) if (eq[k] & RHS) xp |= 1ULL << piv[k];
    int fcols[64], nf = 0;
    for (int c = 0; c < nu; c++) if (!isp[c]) fcols[nf++] = c;
    u64 nv[64];
    for (int f = 0; f < nf; f++) {
        u64 v = 1ULL << fcols[f];
        for (int k = 0; k < npiv; k++) if ((eq[k] >> fcols[f]) & 1) v |= 1ULL << piv[k];
        nv[f] = v;
    }
    if (nf > 36) { fprintf(stderr, "level %d: free dim %d too large\n", j, nf); exit(2); }
    u64 x = xp;
    u64 total = 1ULL << nf;
    for (u64 g = 0; g < total; g++) {
        /* test candidate x */
        int ok = 1;
        u16 Q = JM;
        for (int l = 1; l < M_; l++) {
            u16 t = (u16)((x >> (16 * (l - 1))) & 0xFFFF);
            if (!GL[t]) { ok = 0; break; }
            Q ^= FORM[t];
        }
        if (ok && GL[Q]) {
            for (int l = 1; l < M_; l++) th[j][l] = (u16)((x >> (16 * (l - 1))) & 0xFFFF);
            th[j][0] = I4;
            if (check_conds(j)) {
                nodes[j]++;
                level(j + 1);
                if (maxsol >= 0 && nsol >= maxsol) return;
            }
        }
        if (g + 1 < total) x ^= nv[__builtin_ctzll(g + 1)];
    }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: frames m t v0file classindex [maxsol] [outfile]\n"); return 1; }
    M_ = atoi(argv[1]); T_ = atoi(argv[2]);
    const char *v0file = argv[3]; int ci = atoi(argv[4]);
    if (argc > 5) maxsol = atoll(argv[5]);
    if (argc > 6 && strcmp(argv[6], "-")) fout = fopen(argv[6], "a");
    init_tables();
    I4 = 0; for (int i = 0; i < 4; i++) I4 |= 1 << (5 * i);
    build_conds();
    FILE *f = fopen(v0file, "r");
    if (!f) { perror("v0file"); return 1; }
    char line[512]; int idx = 0;
    while (fgets(line, sizeof line, f)) {
        long sz; unsigned ps[8]; int k = 0;
        char *p = line; char *end;
        sz = strtol(p, &end, 10); p = end;
        for (int l = 1; l < M_; l++) { ps[l] = (unsigned)strtol(p, &end, 10); p = end; k++; }
        if (k != M_ - 1) continue;
        if (ci >= 0 && idx != ci) { idx++; continue; }
        memset(nodes, 0, sizeof nodes); nsol = 0;
        th[0][0] = I4;
        for (int l = 1; l < M_; l++) th[0][l] = (u16)ps[l];
        clock_t t0 = clock();
        int ok = 1;
        u16 Q = JM; for (int l = 1; l < M_; l++) Q ^= FORM[th[0][l]];
        if (!GL[Q]) ok = 0;
        if (ok && !check_conds(0)) ok = 0;
        if (ok) { nodes[0] = 1; level(1); }
        printf("class %d (tuple-orbit %ld): nodes", idx, sz);
        for (int j = 0; j < M_; j++) printf(" %lld", nodes[j]);
        printf("  nsol=%lld  time=%.1fs\n", nsol, (double)(clock() - t0) / CLOCKS_PER_SEC);
        fflush(stdout);
        idx++;
    }
    if (fout) fclose(fout);
    return 0;
}
