/* Completeness test of enum_sub (vsearch8.c): number of 4-dim subspaces inside a universe of U-bit vectors
   must equal the Gaussian binomial [U,4]_2 (U=6: 651, U=7: 11811, U=8: 200787); and with a random candidate set
   the count must equal an independent brute-force count (ordered bases / 20160). */
#define main vsearch_main
#include "vsearch8.c"
#undef main

static long long cnt;
static void count_cb(const int *b, void *ctx) { (void)b; (void)ctx; cnt++; }

int main(void) {
    for (int U = 6; U <= 8; U++) {
        memset(in_c, 0, sizeof in_c);
        int cand[256], nc = 0;
        for (int x = 1; x < (1 << U); x++) { in_c[x] = 1; cand[nc++] = x; }
        cnt = 0;
        enum_sub(cand, nc, count_cb, NULL);
        printf("U=%d: enumerated %lld subspaces\n", U, cnt);
    }
    /* random candidate set in 8-bit universe: compare with ordered-basis brute force */
    srand(12345);
    for (int trial = 0; trial < 3; trial++) {
        memset(in_c, 0, sizeof in_c);
        int cand[256], nc = 0;
        for (int x = 1; x < 256; x++) if (rand() % 100 < 80) { in_c[x] = 1; cand[nc++] = x; }
        cnt = 0;
        enum_sub(cand, nc, count_cb, NULL);
        /* brute force: count ordered bases (x1..x4) with the full span inside the candidate set */
        long long ob = 0;
        for (int a = 1; a < 256; a++) { if (!in_c[a]) continue;
            for (int b = 1; b < 256; b++) { if (b == a || !in_c[b] || !in_c[a ^ b]) continue;
                for (int c = 1; c < 256; c++) {
                    if (c == a || c == b || c == (a ^ b)) continue;
                    if (!in_c[c] || !in_c[c ^ a] || !in_c[c ^ b] || !in_c[c ^ a ^ b]) continue;
                    int sp[8] = {0, a, b, a ^ b, c, c ^ a, c ^ b, c ^ a ^ b};
                    for (int d = 1; d < 256; d++) {
                        int ok = 1;
                        for (int s = 0; s < 8; s++) {
                            if (d == sp[s]) { ok = 0; break; }
                            if (!in_c[d ^ sp[s]]) { ok = 0; break; }
                        }
                        if (ok && in_c[d]) ob++;
                    }
                }
            }
        }
        printf("random trial %d (|cand| = %d): DFS %lld subspaces ; brute-force ordered bases %lld / 20160 = %.3f\n", trial, nc, cnt, ob, ob / 20160.0);
    }
    return 0;
}
