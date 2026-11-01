/* rescheck.c — rigorous bound  Y + (sigma + f) I  >= 0  (p05 verifier).
 *
 * Input on stdin (binary, native endianness):
 *   int32 n, double sigma (exactly the shift used; sigma >= 0)
 *   n(n+1)/2 pairs of doubles (lo, hi): the lower triangle of Y, row by row (i = 0..n-1,
 *   j = 0..i); each pair is an interval containing the exact rational entry.
 * Output: "OK <f>"  or  "FAIL <reason>".
 *
 * Method.  L := floating-point Cholesky factor of mid(Y) + sigma*I (no rigour needed here; any
 * lower-triangular L would do).  F := Y + sigma*I - L L^T is enclosed entrywise by interval
 * arithmetic: every operation is done in round-to-nearest and widened by one ulp each way with
 * nextafter(), which contains the exact result.  f is an upper bound on the Frobenius norm of F,
 * accumulated with upward widening.  Then  Y + (sigma + f) I = L L^T + (F + f I)  is positive
 * semidefinite, because L L^T >= 0 and ||F||_2 <= ||F||_F <= f.
 * Compile: cc -O2 -ffp-contract=off -fno-fast-math -o rescheck rescheck.c -lm
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static inline double dn(double x) { return nextafter(x, -INFINITY); }
static inline double up(double x) { return nextafter(x, INFINITY); }

int main(void) {
    int32_t n;
    double sigma;
    if (fread(&n, sizeof n, 1, stdin) != 1 || n <= 0 || n > 20000) { printf("FAIL bad header\n"); return 2; }
    if (fread(&sigma, sizeof sigma, 1, stdin) != 1 || !(sigma >= 0) || isinf(sigma)) { printf("FAIL bad sigma\n"); return 2; }
    size_t m = (size_t)n * (size_t)(n + 1) / 2;
    double *lo = malloc(m * sizeof(double)), *hi = malloc(m * sizeof(double)), *L = malloc(m * sizeof(double));
    if (!lo || !hi || !L) { printf("FAIL out of memory\n"); return 2; }
    for (size_t t = 0; t < m; t++) {
        double p[2];
        if (fread(p, sizeof(double), 2, stdin) != 2) { printf("FAIL short input\n"); return 2; }
        if (isnan(p[0]) || isnan(p[1]) || isinf(p[0]) || isinf(p[1]) || p[0] > p[1]) { printf("FAIL bad interval\n"); return 1; }
        lo[t] = p[0]; hi[t] = p[1];
    }
#define IX(i, j) ((size_t)(i) * ((i) + 1) / 2 + (j))
    /* 1. floating Cholesky of mid(Y) + sigma I (row-oriented, lower triangle) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double s = 0.5 * (lo[IX(i, j)] + hi[IX(i, j)]) + (i == j ? sigma : 0.0);
            for (int k = 0; k < j; k++) s -= L[IX(i, k)] * L[IX(j, k)];
            if (i == j) {
                if (!(s > 0)) { printf("FAIL cholesky pivot %d not positive\n", i); return 1; }
                L[IX(i, i)] = sqrt(s);
            } else {
                L[IX(i, j)] = s / L[IX(j, j)];
            }
            if (isnan(L[IX(i, j)]) || isinf(L[IX(i, j)])) { printf("FAIL non-finite factor\n"); return 1; }
        }
    }
    /* 2. enclosure of F = Y + sigma I - L L^T and an upward bound on sum F_ij^2 */
    double ss = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double flo = lo[IX(i, j)], fhi = hi[IX(i, j)];
            if (i == j) { flo = dn(flo + sigma); fhi = up(fhi + sigma); }
            for (int k = 0; k <= j; k++) {
                double p = L[IX(i, k)] * L[IX(j, k)];
                double plo = dn(p), phi = up(p);       /* contains the exact product */
                flo = dn(flo - phi);
                fhi = up(fhi - plo);
            }
            double mag = fabs(flo) > fabs(fhi) ? fabs(flo) : fabs(fhi);
            double sq = up(mag * mag);
            if (i != j) sq = up(2.0 * sq);
            ss = up(ss + sq);
            if (isinf(ss) || isnan(ss)) { printf("FAIL residual overflow\n"); return 1; }
        }
    }
    double f = up(sqrt(ss));
    printf("OK %.17g\n", f);
    free(lo); free(hi); free(L);
    return 0;
}
