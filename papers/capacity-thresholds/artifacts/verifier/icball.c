/* icball.c -- rigorous block spectra for p08 (coherent information of permutation-symmetric
 * inputs to n uses of a qubit Pauli channel).  Called by verify.py; see PROBLEM.md "Block formula".
 *
 * Every computed quantity is a ball (centre, radius) with |true - centre| <= radius, kept with
 * round-to-nearest arithmetic plus explicit error terms (u = 2^-53 per rounding, eta = 2^-1074
 * for underflow) and upward-rounded radius arithmetic (up() = next double towards +inf).
 *
 * Path 0 (depolarizing):  block_k = s_k * sum_s c_s E^s Tr^{2k+s}(Y)       (Horner in s)
 * Path 1 (any Pauli):     block_k = s_k * D_d^{-1/2} DP_d( det^k * (Y_hat) ) D_d^{-1/2}
 * Input (stdin, text; doubles in C99 hex):
 *   n r path
 *   path 0: for k = 0..n/2: d+1 coefficient balls c_0..c_d          (mid rad pairs)
 *   path 1: cA cB cC cE balls; det balls D00 D11 D22 D02 D20; invsqrt C(n,a) balls a=0..n;
 *           for k = 0..n/2: invsqrt C(d,m) balls m=0..d
 *   for k = 0..n/2: scale ball s_k
 *   r*(n+1) amplitudes "re im" (exact doubles; row-major: row q, index a = number of ones)
 * Output, per block k:
 *   BLOCK k d
 *   TR S mid rad / TR J mid rad          trace balls of sigma_k and tau_k
 *   EIG S m emb delta  then m eigenvalue centres (descending)    [S = sigma, J = tau]
 *   EIG J ...
 * with |lambda_i(true) - centre_i| <= delta (Weyl pairing of sorted lists).  emb = 1 means the
 * Hermitian matrix was embedded as a real symmetric matrix of twice the size (each eigenvalue
 * listed twice).  Any failure prints "FAIL reason" and exits 2.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define UR 0x1p-53
#define ETA 0x1p-1074

static void fail(const char *msg) {
    printf("FAIL %s\n", msg);
    fflush(stdout);
    exit(2);
}

static void *xmalloc(size_t n) {
    void *p = calloc(1, n ? n : 1);
    if (!p) fail("out of memory");
    return p;
}

/* smallest double >= x + something positive: an upward-rounding step */
static inline double up(double x) {
    if (x == 0.0) return ETA;
    if (x > 0.0 && x < INFINITY) {
        uint64_t b;
        memcpy(&b, &x, 8);
        b += 1;
        memcpy(&x, &b, 8);
        return x;
    }
    return nextafter(x, INFINITY);
}
static inline double dn(double x) { return nextafter(x, -INFINITY); }
static inline double addu(double a, double b) { return up(a + b); }
static inline double mulu(double a, double b) { return up(a * b); } /* a, b >= 0 */

typedef struct { double m, r; } rb;       /* real ball */
typedef struct { double re, im, r; } cb;  /* complex ball */

static inline double abs1(cb a) { return addu(fabs(a.re), fabs(a.im)); } /* >= |centre| */

static inline cb cb_add(cb a, cb b) {
    cb c;
    c.re = a.re + b.re;
    c.im = a.im + b.im;
    c.r = addu(addu(a.r, b.r), mulu(UR, addu(fabs(c.re), fabs(c.im))));
    return c;
}

static inline cb cb_smul(rb s, cb a) { /* real ball times complex ball */
    cb c;
    c.re = s.m * a.re;
    c.im = s.m * a.im;
    double err = addu(mulu(UR, addu(fabs(c.re), fabs(c.im))), ETA);
    c.r = addu(addu(addu(mulu(fabs(s.m), a.r), mulu(s.r, abs1(a))), mulu(s.r, a.r)), err);
    return c;
}

static inline cb cb_mul(cb a, cb b) {
    cb c;
    double p1 = a.re * b.re, p2 = a.im * b.im, p3 = a.re * b.im, p4 = a.im * b.re;
    c.re = p1 - p2;
    c.im = p3 + p4;
    double a1 = abs1(a), b1 = abs1(b);
    double err = addu(mulu(3 * UR, mulu(a1, b1)), 4 * ETA);
    c.r = addu(addu(addu(mulu(a1, b.r), mulu(b1, a.r)), mulu(a.r, b.r)), err);
    return c;
}

static inline cb cb_conj(cb a) { a.im = -a.im; return a; }
static const cb CB0 = {0.0, 0.0, 0.0};

/* sqrt(x)/N for exact integers x < 2^53, N > 0: correctly rounded sqrt then division */
static inline rb wsq(double x, double N) {
    rb w;
    w.m = sqrt(x) / N;
    w.r = mulu(3 * UR, w.m);
    return w;
}

/* ---------------- input ---------------- */

static double rd(void) {
    char buf[128];
    if (scanf("%127s", buf) != 1) fail("short input");
    char *end;
    double v = strtod(buf, &end);
    if (*end) fail("bad number");
    if (!isfinite(v)) fail("non-finite input");
    return v;
}
static rb rdb(void) {
    rb b;
    b.m = rd();
    b.r = rd();
    if (!(b.r >= 0)) fail("negative radius");
    return b;
}

/* ---------------- path 0: depolarizing ---------------- */

/* T (size N+1) -> out (size N): one-qubit partial trace in the Dicke basis */
static void tr1(const cb *T, int N, cb *out) {
    int S = N + 1;
    for (int c = 0; c < N; c++)
        for (int d = 0; d < N; d++) {
            rb w0 = wsq((double)(N - c) * (double)(N - d), N);
            rb w1 = wsq((double)(c + 1) * (double)(d + 1), N);
            out[c * N + d] = cb_add(cb_smul(w0, T[c * S + d]), cb_smul(w1, T[(c + 1) * S + d + 1]));
        }
}

/* Z (size N) -> out (size N+1): E = Tr_1^dagger, P_sym (Z (x) I) P_sym */
static void emb1(const cb *Z, int N, cb *out) {
    int S = N + 1;
    for (int a = 0; a <= N; a++)
        for (int b = 0; b <= N; b++) {
            cb acc = CB0;
            if (a < N && b < N) acc = cb_smul(wsq((double)(N - a) * (double)(N - b), N), Z[a * N + b]);
            if (a >= 1 && b >= 1)
                acc = cb_add(acc, cb_smul(wsq((double)a * (double)b, N), Z[(a - 1) * N + b - 1]));
            out[a * S + b] = acc;
        }
}

/* ---------------- eigenvalue certification ---------------- */

/* JAMA-style Householder tridiagonalisation + implicit QL (hint only: soundness never depends on
   it, a bad decomposition just produces a large delta). V is n x n row-major, in: A, out: vecs
   in columns; d: eigenvalues.  returns 0 on success. */
static int eig_sym(int n, double *V, double *d, double *e) {
#define Vx(i, j) V[(size_t)(i) * n + (j)]
    for (int j = 0; j < n; j++) d[j] = Vx(n - 1, j);
    for (int i = n - 1; i > 0; i--) {
        double scale = 0.0, h = 0.0;
        for (int k = 0; k < i; k++) scale += fabs(d[k]);
        if (scale == 0.0) {
            e[i] = d[i - 1];
            for (int j = 0; j < i; j++) {
                d[j] = Vx(i - 1, j);
                Vx(i, j) = 0.0;
                Vx(j, i) = 0.0;
            }
        } else {
            for (int k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            double f = d[i - 1], g = sqrt(h);
            if (f > 0) g = -g;
            e[i] = scale * g;
            h = h - f * g;
            d[i - 1] = f - g;
            for (int j = 0; j < i; j++) e[j] = 0.0;
            for (int j = 0; j < i; j++) {
                f = d[j];
                Vx(j, i) = f;
                g = e[j] + Vx(j, j) * f;
                for (int k = j + 1; k <= i - 1; k++) {
                    g += Vx(k, j) * d[k];
                    e[k] += Vx(k, j) * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (int j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            double hh = f / (h + h);
            for (int j = 0; j < i; j++) e[j] -= hh * d[j];
            for (int j = 0; j < i; j++) {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i - 1; k++) Vx(k, j) -= (f * e[k] + g * d[k]);
                d[j] = Vx(i - 1, j);
                Vx(i, j) = 0.0;
            }
        }
        d[i] = h;
    }
    for (int i = 0; i < n - 1; i++) {
        Vx(n - 1, i) = Vx(i, i);
        Vx(i, i) = 1.0;
        double h = d[i + 1];
        if (h != 0.0) {
            for (int k = 0; k <= i; k++) d[k] = Vx(k, i + 1) / h;
            for (int j = 0; j <= i; j++) {
                double g = 0.0;
                for (int k = 0; k <= i; k++) g += Vx(k, i + 1) * Vx(k, j);
                for (int k = 0; k <= i; k++) Vx(k, j) -= g * d[k];
            }
        }
        for (int k = 0; k <= i; k++) Vx(k, i + 1) = 0.0;
    }
    for (int j = 0; j < n; j++) {
        d[j] = Vx(n - 1, j);
        Vx(n - 1, j) = 0.0;
    }
    Vx(n - 1, n - 1) = 1.0;
    e[0] = 0.0;
    /* tql2 */
    for (int i = 1; i < n; i++) e[i - 1] = e[i];
    e[n - 1] = 0.0;
    double f = 0.0, tst1 = 0.0, eps = 0x1p-52;
    for (int l = 0; l < n; l++) {
        tst1 = fmax(tst1, fabs(d[l]) + fabs(e[l]));
        int m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps * tst1) break;
            m++;
        }
        if (m == n) m = n - 1;
        if (m > l) {
            int iter = 0;
            do {
                if (++iter > 200) return 1;
                double g = d[l];
                double p = (d[l + 1] - g) / (2.0 * e[l]);
                double r = hypot(p, 1.0);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                double dl1 = d[l + 1];
                double h = g - d[l];
                for (int i = l + 2; i < n; i++) d[i] -= h;
                f = f + h;
                p = d[m];
                double c = 1.0, c2 = c, c3 = c, el1 = e[l + 1], s = 0.0, s2 = 0.0;
                for (int i = m - 1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for (int k = 0; k < n; k++) {
                        h = Vx(k, i + 1);
                        Vx(k, i + 1) = s * Vx(k, i) + c * h;
                        Vx(k, i) = c * Vx(k, i) - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
            } while (fabs(e[l]) > eps * tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }
    for (int i = 0; i < n; i++)
        if (!isfinite(d[i])) return 1;
    return 0;
#undef Vx
}

/* Dot2 (Ogita-Rump-Oishi 2005): x.y + extra, computed as if in twice the working precision.
   |res - exact| <= u|exact| + gamma_{len}^2 sum|x_i y_i| (no underflow; underflow is absorbed by
   the explicit eta terms added by the caller). */
static inline double dot2(const double *x, const double *y, int n, double xe, double ye) {
    double p = xe * ye, s = fma(xe, ye, -p);
    for (int i = 0; i < n; i++) {
        double h = x[i] * y[i];
        double r = fma(x[i], y[i], -h);
        double t = p + h;
        double z = t - p;
        double q = (p - (t - z)) + (h - z);
        p = t;
        s += q + r;
    }
    return p + s;
}

static int cmp_desc(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) - (x > y);
}

/* certify the spectrum of the Hermitian ball matrix M (m x m, upper triangle used) */
static void certify(const char *tag, int m, const cb *M) {
    /* trace ball, radius Frobenius norm, and whether the centre is real */
    double trm = 0.0, trr = 0.0, rad2 = 0.0;
    int cplx = 0;
    for (int i = 0; i < m; i++) {
        double t = trm + M[(size_t)i * m + i].re;
        trr = addu(addu(trr, M[(size_t)i * m + i].r), mulu(UR, fabs(t)));
        trm = t;
        for (int j = i; j < m; j++) {
            cb x = M[(size_t)i * m + j];
            double q = mulu(x.r, x.r);
            rad2 = addu(rad2, j == i ? q : mulu(2.0, q));
            if (j > i && x.im != 0.0) cplx = 1;
        }
    }
    printf("TR %s %a %a\n", tag, trm, trr);
    double radF = up(sqrt(rad2));
    int n = cplx ? 2 * m : m;
    double *A = xmalloc(sizeof(double) * (size_t)n * n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) {
            cb x = i <= j ? M[(size_t)i * m + j] : cb_conj(M[(size_t)j * m + i]);
            double re = x.re, im = (i == j) ? 0.0 : x.im;
            A[(size_t)i * n + j] = re;
            if (cplx) {
                A[(size_t)(i + m) * n + (j + m)] = re;
                A[(size_t)i * n + (j + m)] = -im;
                A[(size_t)(i + m) * n + j] = im;
            }
        }
    double *V = xmalloc(sizeof(double) * (size_t)n * n);
    memcpy(V, A, sizeof(double) * (size_t)n * n);
    double *lam = xmalloc(sizeof(double) * n), *e = xmalloc(sizeof(double) * n);
    if (eig_sym(n, V, lam, e)) fail("eigensolver did not converge");
    /* QT: eigenvectors as rows */
    double *QT = xmalloc(sizeof(double) * (size_t)n * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) QT[(size_t)j * n + i] = V[(size_t)i * n + j];
    free(V);
    /* norms (upward) */
    double nA2 = 0.0, nQ2 = 0.0, lmax = 0.0;
    for (size_t t = 0; t < (size_t)n * n; t++) {
        nA2 = addu(nA2, mulu(fabs(A[t]), fabs(A[t])));
        nQ2 = addu(nQ2, mulu(fabs(QT[t]), fabs(QT[t])));
    }
    for (int i = 0; i < n; i++) lmax = fmax(lmax, fabs(lam[i]));
    double nA = up(sqrt(nA2)), nQ = up(sqrt(nQ2));
    double k1 = (double)n + 1.0;
    double gam = up(up(k1 * UR) / dn(1.0 - k1 * UR)); /* gamma_{n+1} */
    double g2 = mulu(gam, gam);
    /* residual W = A Q - Q Lambda and orthogonality E = Q^T Q - I, entrywise by Dot2 */
    double w2 = 0.0, e2 = 0.0;
    for (int j = 0; j < n; j++) {
        const double *qj = QT + (size_t)j * n;
        for (int i = 0; i < n; i++) {
            double w = dot2(A + (size_t)i * n, qj, n, -lam[j], qj[i]);
            w2 = addu(w2, mulu(fabs(w), fabs(w)));
        }
        for (int i = 0; i < n; i++) {
            double v = dot2(QT + (size_t)i * n, qj, n, i == j ? -1.0 : 0.0, 1.0);
            e2 = addu(e2, mulu(fabs(v), fabs(v)));
        }
    }
    /* ||W_true||_F <= (||W_comp||_F + g2 ||S||_F)/(1-u),  ||S||_F <= (||A||_F + lmax) ||Q||_F */
    double uf = dn(1.0 - UR);
    double tiny = mulu(mulu(4.0 * ETA, k1), k1);
    double sW = mulu(addu(nA, lmax), nQ);
    double nW = up(addu(addu(up(sqrt(w2)), mulu(g2, sW)), tiny) / uf);
    double sE = addu(mulu(nQ, nQ), up(sqrt((double)n)));
    double eta = up(addu(addu(up(sqrt(e2)), mulu(g2, sE)), tiny) / uf);
    if (!(eta < 0.5)) fail("eigenvectors not orthonormal enough");
    /* Weyl via polar factor: |lam_i(A) - lam_i| <= (||W|| + 2 eta ||Lam||)/(1 - eta) */
    double delta = up(addu(nW, mulu(2.0 * eta, lmax)) / dn(1.0 - eta));
    delta = addu(delta, radF); /* ||M_true - M_centre||_2 <= ||rad||_F, embedding keeps the norm */
    if (!isfinite(delta)) fail("non-finite delta");
    qsort(lam, n, sizeof(double), cmp_desc);
    printf("EIG %s %d %d %a\n", tag, n, cplx, delta);
    for (int i = 0; i < n; i++) printf("%a\n", lam[i]);
    free(A);
    free(QT);
    free(lam);
    free(e);
}

/* ---------------- main ---------------- */

int main(void) {
    int n, r, path;
    if (scanf("%d %d %d", &n, &r, &path) != 3) fail("bad header");
    if (n < 1 || n > 1000 || r < 1 || r > 8 || (path != 0 && path != 1)) fail("bad header values");
    int K = n / 2;
    rb **coef = NULL;
    rb cA = {0, 0}, cB = {0, 0}, cC = {0, 0}, cE = {0, 0}, D[5];
    rb *isn = NULL, **isd = NULL;
    if (path == 0) {
        coef = xmalloc(sizeof(rb *) * (K + 1));
        for (int k = 0; k <= K; k++) {
            int d = n - 2 * k;
            coef[k] = xmalloc(sizeof(rb) * (d + 1));
            for (int s = 0; s <= d; s++) coef[k][s] = rdb();
        }
    } else {
        cA = rdb(); cB = rdb(); cC = rdb(); cE = rdb();
        for (int i = 0; i < 5; i++) D[i] = rdb();
        isn = xmalloc(sizeof(rb) * (n + 1));
        for (int a = 0; a <= n; a++) isn[a] = rdb();
        isd = xmalloc(sizeof(rb *) * (K + 1));
        for (int k = 0; k <= K; k++) {
            int d = n - 2 * k;
            isd[k] = xmalloc(sizeof(rb) * (d + 1));
            for (int m = 0; m <= d; m++) isd[k][m] = rdb();
        }
    }
    rb *scale = xmalloc(sizeof(rb) * (K + 1));
    for (int k = 0; k <= K; k++) scale[k] = rdb();
    int S0 = n + 1;
    cb *amp = xmalloc(sizeof(cb) * r * S0);
    for (int q = 0; q < r; q++)
        for (int a = 0; a <= n; a++) {
            amp[q * S0 + a].re = rd();
            amp[q * S0 + a].im = rd();
            amp[q * S0 + a].r = 0.0;
        }
    /* blocks[k][q1*r+q2] : (d+1)^2 balls, only q1 <= q2 computed */
    cb ***blk = xmalloc(sizeof(cb **) * (K + 1));
    for (int k = 0; k <= K; k++) blk[k] = xmalloc(sizeof(cb *) * r * r);
    cb *Y = xmalloc(sizeof(cb) * S0 * S0);
    for (int q1 = 0; q1 < r; q1++)
        for (int q2 = q1; q2 < r; q2++) {
            for (int a = 0; a <= n; a++)
                for (int b = 0; b <= n; b++)
                    Y[a * S0 + b] = cb_mul(amp[q1 * S0 + a], cb_conj(amp[q2 * S0 + b]));
            if (path == 0) {
                /* chain T_t = Tr^t Y, t = 0..n, sizes n-t+1 */
                cb **T = xmalloc(sizeof(cb *) * (n + 1));
                T[0] = xmalloc(sizeof(cb) * S0 * S0);
                memcpy(T[0], Y, sizeof(cb) * S0 * S0);
                for (int t = 1; t <= n; t++) {
                    int N = n - t + 1;
                    T[t] = xmalloc(sizeof(cb) * (size_t)N * N);
                    tr1(T[t - 1], N, T[t]);
                }
                for (int k = 0; k <= K; k++) {
                    int d = n - 2 * k;
                    cb *R = xmalloc(sizeof(cb));
                    R[0] = cb_smul(coef[k][d], T[n][0]);
                    for (int s = d - 1; s >= 0; s--) {
                        int L = d - s + 1; /* size of T_{2k+s} */
                        cb *R2 = xmalloc(sizeof(cb) * (size_t)L * L);
                        emb1(R, L - 1, R2);
                        const cb *Tt = T[2 * k + s];
                        for (int i = 0; i < L * L; i++) R2[i] = cb_add(R2[i], cb_smul(coef[k][s], Tt[i]));
                        free(R);
                        R = R2;
                    }
                    for (int i = 0; i < (d + 1) * (d + 1); i++) R[i] = cb_smul(scale[k], R[i]);
                    blk[k][q1 * r + q2] = R;
                }
                for (int t = 0; t <= n; t++) free(T[t]);
                free(T);
            } else {
                /* monomial coefficients, then det-contraction chain and the DP per block */
                int sz = S0;
                cb *Yk = xmalloc(sizeof(cb) * S0 * S0);
                for (int a = 0; a <= n; a++)
                    for (int b = 0; b <= n; b++) {
                        cb x = cb_smul(isn[a], Y[a * S0 + b]);
                        Yk[a * S0 + b] = cb_smul(isn[b], x);
                    }
                static const int dal[5] = {0, 1, 2, 0, 2}, dbe[5] = {0, 1, 2, 2, 0};
                for (int k = 0; k <= K; k++) {
                    if (k > 0) {
                        int m = sz - 2;
                        cb *nx = xmalloc(sizeof(cb) * (size_t)m * m);
                        for (int a = 0; a < m; a++)
                            for (int b = 0; b < m; b++) {
                                cb acc = CB0;
                                for (int t = 0; t < 5; t++)
                                    acc = cb_add(acc, cb_smul(D[t], Yk[(a + dal[t]) * sz + b + dbe[t]]));
                                nx[a * m + b] = acc;
                            }
                        free(Yk);
                        Yk = nx;
                        sz = m;
                    }
                    int d = n - 2 * k; /* sz == d+1 */
                    /* G[(al*(s+1)+be)*L*L + mu*L + nu], level s, L = d-s+1; two ping-pong
                       buffers of the largest level size (sizes first grow, then shrink) */
                    int s = d, L = 1;
                    size_t mx = 0;
                    for (int t = 0; t <= d; t++) {
                        size_t z = (size_t)(t + 1) * (t + 1) * (d - t + 1) * (d - t + 1);
                        if (z > mx) mx = z;
                    }
                    cb *G = xmalloc(sizeof(cb) * mx), *Gb = xmalloc(sizeof(cb) * mx);
                    for (int al = 0; al <= d; al++)
                        for (int be = 0; be <= d; be++) G[al * (d + 1) + be] = Yk[al * sz + be];
                    /* kernels: (da, db) -> two (dm, dm') shifts with weights */
                    const int kda[4] = {0, 1, 0, 1}, kdb[4] = {0, 1, 1, 0};
                    const int kdm[4][2][2] = {{{0, 0}, {1, 1}}, {{0, 0}, {1, 1}}, {{0, 1}, {1, 0}}, {{1, 0}, {0, 1}}};
                    const rb kw[4][2] = {{cA, cB}, {cB, cA}, {cC, cE}, {cC, cE}};
                    while (s > 0) {
                        int ns = s - 1, nL = L + 1, w1 = s + 1, w0 = ns + 1;
                        cb *G2 = Gb;
                        memset(G2, 0, sizeof(cb) * (size_t)w0 * w0 * nL * nL);
                        for (int al = 0; al <= ns; al++)
                            for (int be = 0; be <= ns; be++) {
                                cb *out = G2 + ((size_t)al * w0 + be) * nL * nL;
                                for (int t = 0; t < 4; t++) {
                                    const cb *X = G + ((size_t)(al + kda[t]) * w1 + be + kdb[t]) * L * L;
                                    for (int v = 0; v < 2; v++) {
                                        rb w = kw[t][v];
                                        if (w.m == 0.0 && w.r == 0.0) continue;
                                        int dm = kdm[t][v][0], dm2 = kdm[t][v][1];
                                        for (int mu = 0; mu < L; mu++)
                                            for (int nu = 0; nu < L; nu++) {
                                                cb *o = out + (mu + dm) * nL + nu + dm2;
                                                *o = cb_add(*o, cb_smul(w, X[mu * L + nu]));
                                            }
                                    }
                                }
                            }
                        Gb = G;
                        G = G2;
                        s = ns;
                        L = nL;
                    }
                    free(Gb);
                    /* G is now level 0: single array L x L, L = d+1 */
                    cb *Bk = xmalloc(sizeof(cb) * (size_t)L * L);
                    for (int mu = 0; mu <= d; mu++)
                        for (int nu = 0; nu <= d; nu++) {
                            cb x = cb_smul(isd[k][mu], G[mu * L + nu]);
                            x = cb_smul(isd[k][nu], x);
                            Bk[mu * L + nu] = cb_smul(scale[k], x);
                        }
                    free(G);
                    blk[k][q1 * r + q2] = Bk;
                }
                free(Yk);
            }
        }
    /* assemble and certify each block */
    for (int k = 0; k <= K; k++) {
        int d = n - 2 * k, D1 = d + 1, m = r * D1;
        printf("BLOCK %d %d\n", k, d);
        cb *Sg = xmalloc(sizeof(cb) * D1 * D1);
        for (int q = 0; q < r; q++) {
            const cb *B = blk[k][q * r + q];
            for (int i = 0; i < D1 * D1; i++) Sg[i] = q ? cb_add(Sg[i], B[i]) : B[i];
        }
        certify("S", D1, Sg);
        free(Sg);
        cb *J = xmalloc(sizeof(cb) * (size_t)m * m);
        for (int q1 = 0; q1 < r; q1++)
            for (int q2 = q1; q2 < r; q2++) {
                const cb *B = blk[k][q1 * r + q2];
                for (int a = 0; a < D1; a++)
                    for (int b = 0; b < D1; b++) {
                        J[(size_t)(q1 * D1 + a) * m + q2 * D1 + b] = B[a * D1 + b];
                        if (q2 > q1) J[(size_t)(q2 * D1 + b) * m + q1 * D1 + a] = cb_conj(B[a * D1 + b]);
                    }
            }
        certify("J", m, J);
        free(J);
        for (int q1 = 0; q1 < r; q1++)
            for (int q2 = q1; q2 < r; q2++) free(blk[k][q1 * r + q2]);
        fflush(stdout);
    }
    printf("END\n");
    return 0;
}
