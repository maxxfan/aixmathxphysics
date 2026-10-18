/* p04 verifier: exact integer accumulation of a local-model certificate.
 *
 *   accum IN OUT
 *
 * IN (little-endian, written by verify.py):
 *   int64 mA, mB, K, G, transpose
 *   K records: int64 W, then ceil(mA/64) uint64 words of Alice's signs, ceil(mB/64) words of Bob's
 *              (bit i set = +1, clear = -1; word w holds indices 64w .. 64w+63, bit 0 first)
 *   G records: int32 piA[mA], int8 sgA[mA], int32 piB[mB], int8 sgB[mB]
 * OUT:
 *   int64 N0[mA*mB]  (row-major)  N0[i][j] = sum_k W_k sA_k(i) sB_k(j)
 *   int64 S[mA*mB]                S[i][j]  = sum_h sgA_h(i) sgB_h(j) N0[piA_h(i)][piB_h(j)]
 *                                 (then S + S^T if transpose)
 * The caller guarantees sum_k W_k <= 2^50 and G * (1 + transpose) <= 4096, so every partial sum is
 * bounded by 2^62 in absolute value and int64 arithmetic is exact.  Exit code != 0 on any error.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void die(const char *m) { fprintf(stderr, "accum: %s\n", m); exit(2); }

static void rd(void *p, size_t sz, size_t n, FILE *f) {
    if (fread(p, sz, n, f) != n) die("short read");
}

int main(int argc, char **argv) {
    if (argc != 3) die("usage: accum IN OUT");
    FILE *f = fopen(argv[1], "rb");
    if (!f) die("cannot open input");
    int64_t hdr[5];
    rd(hdr, 8, 5, f);
    int64_t mA = hdr[0], mB = hdr[1], K = hdr[2], G = hdr[3], tr = hdr[4];
    if (mA < 1 || mB < 1 || mA > 4096 || mB > 4096 || K < 1 || K > 20000000 || G < 1 || G > 4096 ||
        (tr != 0 && tr != 1) || (tr && mA != mB))
        die("bad header");
    int64_t wa = (mA + 63) / 64, wb = (mB + 63) / 64;
    int64_t *N0 = calloc((size_t)(mA * mB), 8);
    int64_t *tw = malloc((size_t)mB * 8);
    uint64_t *ba = malloc((size_t)wa * 8), *bb = malloc((size_t)wb * 8);
    if (!N0 || !tw || !ba || !bb) die("out of memory");
    for (int64_t k = 0; k < K; k++) {
        int64_t W;
        rd(&W, 8, 1, f);
        rd(ba, 8, (size_t)wa, f);
        rd(bb, 8, (size_t)wb, f);
        if (W <= 0) die("non-positive weight");
        for (int64_t j = 0; j < mB; j++) tw[j] = ((bb[j >> 6] >> (j & 63)) & 1) ? W : -W;
        for (int64_t i = 0; i < mA; i++) {
            int64_t *row = N0 + i * mB;
            if ((ba[i >> 6] >> (i & 63)) & 1)
                for (int64_t j = 0; j < mB; j++) row[j] += tw[j];
            else
                for (int64_t j = 0; j < mB; j++) row[j] -= tw[j];
        }
    }
    int64_t *S = calloc((size_t)(mA * mB), 8);
    int32_t *pa = malloc((size_t)mA * 4), *pb = malloc((size_t)mB * 4);
    int8_t *sa = malloc((size_t)mA), *sb = malloc((size_t)mB);
    if (!S || !pa || !pb || !sa || !sb) die("out of memory");
    for (int64_t h = 0; h < G; h++) {
        rd(pa, 4, (size_t)mA, f);
        rd(sa, 1, (size_t)mA, f);
        rd(pb, 4, (size_t)mB, f);
        rd(sb, 1, (size_t)mB, f);
        for (int64_t i = 0; i < mA; i++)
            if (pa[i] < 0 || pa[i] >= mA || (sa[i] != 1 && sa[i] != -1)) die("bad Alice transform");
        for (int64_t j = 0; j < mB; j++)
            if (pb[j] < 0 || pb[j] >= mB || (sb[j] != 1 && sb[j] != -1)) die("bad Bob transform");
        for (int64_t i = 0; i < mA; i++) {
            const int64_t *src = N0 + (int64_t)pa[i] * mB;
            int64_t *dst = S + i * mB;
            for (int64_t j = 0; j < mB; j++) {
                int64_t v = src[pb[j]];
                dst[j] += (sa[i] * sb[j] > 0) ? v : -v;
            }
        }
    }
    {
        char extra;
        if (fread(&extra, 1, 1, f) != 0) die("trailing input");
    }
    fclose(f);
    if (tr) {
        for (int64_t i = 0; i < mA; i++)
            for (int64_t j = i + 1; j < mA; j++) {
                int64_t s = S[i * mA + j] + S[j * mA + i];
                S[i * mA + j] = s;
                S[j * mA + i] = s;
            }
        for (int64_t i = 0; i < mA; i++) S[i * mA + i] *= 2;
    }
    FILE *o = fopen(argv[2], "wb");
    if (!o) die("cannot open output");
    if (fwrite(N0, 8, (size_t)(mA * mB), o) != (size_t)(mA * mB)) die("short write");
    if (fwrite(S, 8, (size_t)(mA * mB), o) != (size_t)(mA * mB)) die("short write");
    if (fclose(o) != 0) die("close failed");
    return 0;
}
