/* phasecheck.c -- independent functional check of a T-count certificate for GF(2^n) multiplication.
 *
 * The multiplier |a,b,c> -> |a,b,c + a*b> equals H_c D H_c with D = diag(w^{4 f(x)}),
 * w = exp(i pi/4), f(a,b,c) = c . (a*b mod P).  Factors u_1..u_m give the T-gate phase
 * polynomial sum_r l_r(x), l_r(x) = (u_r . x mod 2) as an integer.  D equals
 * (prod_r T-gadget(u_r)) times a diagonal Clifford iff
 *     G(x) = 4 f(x) - sum_r l_r(x)  (mod 8)
 * is even for all x and h = G/2 (mod 4) has a multilinear expansion over Z_4 with no terms of
 * degree >= 3 and only even coefficients in degree 2 (those are S and CZ gates).
 * sum_r l_r(x) = (m - W(x)) / 2, with W the Walsh-Hadamard transform of the factor multiset.
 *
 * Written independently of verify.py (own field arithmetic, no shared code).
 * input:  n P m  then m factor ints (bit i = variable i; a_i = i, b_j = n+j, c_k = 2n+k)
 * output: one JSON object {"ok": true|false, "reason": "...", "N": .., "seconds": ..}
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static uint32_t gmul(uint32_t a, uint32_t b, uint32_t P, int n) {
    uint64_t r = 0;
    for (int i = 0; i < n; i++)
        if (b >> i & 1) r ^= (uint64_t)a << i;
    for (int k = 2 * n - 2; k >= n; k--)
        if (r >> k & 1) r ^= (uint64_t)P << (k - n);
    return (uint32_t)r;
}

int main(void) {
    int n, m;
    unsigned P;
    if (scanf("%d %u %d", &n, &P, &m) != 3 || n < 1 || n > 9) {
        printf("{\"ok\": false, \"reason\": \"bad header or n > 9\"}\n");
        return 2;
    }
    double t0 = now();
    int N = 3 * n;
    uint64_t size = 1ULL << N, mask = (1ULL << n) - 1;
    int32_t *W = calloc(size, sizeof(int32_t));
    uint8_t *h = malloc(size);
    uint32_t *prod = malloc(sizeof(uint32_t) << (2 * n));
    if (!W || !h || !prod) { printf("{\"ok\": false, \"reason\": \"out of memory\"}\n"); return 2; }
    for (int i = 0; i < m; i++) {
        unsigned long long u;
        if (scanf("%llu", &u) != 1 || u == 0 || u >= size) {
            printf("{\"ok\": false, \"reason\": \"bad factor %d\"}\n", i);
            return 2;
        }
        W[u] += 1;
    }
    /* Walsh-Hadamard: W(x) = sum_v F[v] (-1)^{v.x} */
    for (int b = 0; b < N; b++) {
        uint64_t s = 1ULL << b;
        for (uint64_t x = 0; x < size; x++)
            if (x & s) {
                int32_t p = W[x ^ s], q = W[x];
                W[x ^ s] = p + q;
                W[x] = p - q;
            }
    }
    for (uint32_t a = 0; a <= mask; a++)
        for (uint32_t bb = 0; bb <= mask; bb++) prod[(a << n) | bb] = gmul(a, bb, P, n);
    for (uint64_t x = 0; x < size; x++) {
        uint32_t a = x & mask, bb = (x >> n) & mask, c = (uint32_t)(x >> (2 * n));
        int f = __builtin_popcount(c & prod[(a << n) | bb]) & 1;
        int sumL = (m - W[x]) / 2;
        int G = ((4 * f - sumL) % 8 + 8) % 8;
        if (G & 1) {
            printf("{\"ok\": false, \"reason\": \"phase difference is odd at x=%llu (not Clifford)\", \"N\": %d}\n",
                   (unsigned long long)x, N);
            return 0;
        }
        h[x] = (uint8_t)((G / 2) & 3);
    }
    /* Moebius transform over Z_4: h_S = sum_{T subset S} (-1)^{|S|-|T|} h(1_T) */
    for (int b = 0; b < N; b++) {
        uint64_t s = 1ULL << b;
        for (uint64_t x = 0; x < size; x++)
            if (x & s) h[x] = (uint8_t)((h[x] - h[x ^ s]) & 3);
    }
    for (uint64_t x = 0; x < size; x++) {
        int k = __builtin_popcountll(x);
        if ((k >= 3 && h[x]) || (k == 2 && (h[x] & 1))) {
            printf("{\"ok\": false, \"reason\": \"non-Clifford term of degree %d (coefficient %d) on variables %llu\", \"N\": %d}\n",
                   k, h[x], (unsigned long long)x, N);
            return 0;
        }
    }
    printf("{\"ok\": true, \"reason\": \"\", \"N\": %d, \"seconds\": %.3f}\n", N, now() - t0);
    return 0;
}
