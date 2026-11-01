/* fpenv.c — check the floating-point assumptions of verifier/rescheck.c on this machine and compiler.
 *
 * rescheck.c encloses every rounded operation by widening its result one ulp each way with nextafter().
 * That is valid under IEEE 754 binary64 arithmetic with round-to-nearest, evaluation in binary64
 * (FLT_EVAL_METHOD == 0), gradual underflow (no flush-to-zero), no fused multiply-add contraction and a
 * correctly rounded sqrt.  This program tests each of these properties.  Build it with the flags that
 * verify.py uses for rescheck.c:   cc -O2 -ffp-contract=off -fno-fast-math -o fpenv fpenv.c -lm
 * The replay also builds it with -ffp-contract=fast as a negative control: on hardware with FMA (arm64,
 * x86-64 with -mfma) that build reports "contraction FUSED".
 */
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

int main(void) {
    int bad = 0;
    printf("FLT_EVAL_METHOD %d\n", (int)FLT_EVAL_METHOD);
    if (FLT_EVAL_METHOD != 0) bad = 1;
    int rn = fegetround() == FE_TONEAREST;
    printf("rounding %s\n", rn ? "nearest" : "OTHER");
    if (!rn) bad = 1;
    volatile double m = DBL_MIN;
    volatile double h = m / 2;
    int sub = h > 0 && h * 2 == m && nextafter(0.0, 1.0) == 0x1p-1074;
    printf("subnormals %s\n", sub ? "kept" : "FLUSHED");
    if (!sub) bad = 1;
    /* a*b = 1 - 2^-54 exactly, a tie that rounds to 1.0: r = 0 if the product is rounded, -2^-54 if fused */
    volatile double a = 1 + 0x1p-27, b = 1 - 0x1p-27, c = 1.0;
    double r = a * b - c;
    printf("contraction %s\n", r == 0 ? "none" : "FUSED");
    if (r != 0) bad = 1;
    int nx = nextafter(1.0, 2.0) == 1 + 0x1p-52 && nextafter(1.0, 0.0) == 1 - 0x1p-53 &&
             nextafter(0.0, -1.0) == -0x1p-1074 && nextafter(DBL_MAX, INFINITY) == INFINITY;
    printf("nextafter %s\n", nx ? "ok" : "WRONG");
    if (!nx) bad = 1;
    volatile double two = 2.0, three = 3.0;
    int sq = sqrt(two) == 0x1.6a09e667f3bcdp+0 && sqrt(three) == 0x1.bb67ae8584caap+0;
    printf("sqrt %s\n", sq ? "ok" : "WRONG");
    if (!sq) bad = 1;
    volatile double big = DBL_MAX;
    double ov = big * 2;
    printf("overflow %s\n", isinf(ov) ? "inf" : "WRONG");
    if (!isinf(ov)) bad = 1;
    printf("%s\n", bad ? "FPENV FAIL" : "FPENV OK");
    return bad;
}
