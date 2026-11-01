/* DFS over symmetric 8x8 F_2 matrices beta with good[w][lambda(w)] for all w != 0.  Input file: 65536 bytes good[w*256+lam]. Output: "nodes sols beta0..beta7(first)". */
#include <stdio.h>
#include <stdlib.h>
static unsigned char good[256][256];
static int beta[8];
static long long nodes = 0, nsol = 0, limit = 1;
static int first[8], have_first = 0;
static void rec(int k) {
    if (nsol >= limit) return;
    int low = 0;
    for (int j = 0; j < k; j++) low |= ((beta[j] >> k) & 1) << j;
    for (int fr = 0; fr < (1 << (8 - k)); fr++) {
        int row = low | (fr << k);
        beta[k] = row;
        int ok = 1;
        for (int s = 0; s < (1 << k) && ok; s++) {
            int w = (1 << k) | s, lam = 0, ww = w;
            for (int i = 0; ww; i++, ww >>= 1) if (ww & 1) lam ^= beta[i];
            if (!good[w][lam]) ok = 0;
        }
        nodes++;
        if (!ok) continue;
        if (k == 7) {
            nsol++;
            if (!have_first) { for (int i = 0; i < 8; i++) first[i] = beta[i]; have_first = 1; }
            if (nsol >= limit) return;
        } else rec(k + 1);
    }
    beta[k] = 0;
}
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb");
    if (!f || fread(good, 1, 65536, f) != 65536) return 2;
    fclose(f);
    if (argc > 2) limit = atoll(argv[2]);
    rec(0);
    printf("%lld %lld", nodes, nsol);
    if (have_first) for (int i = 0; i < 8; i++) printf(" %d", first[i]);
    printf("\n");
    return 0;
}
