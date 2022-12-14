#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdarg.h>
#include "rng.h"
#include "dme.h"

fq_elem fq_add(fq_elem a, fq_elem b)
{
    return a ^ b;
}

#ifdef REFERENCE
fq_elem fq_mul(fq_elem a, fq_elem b)
{
    unsigned int i;
    fq_elem c;
    c = fq_0;
    for (i=0; i<64; i++)
    {
         c = (c << 1) ^ ((-(c >> 63)) & fq_t64) ^ ((-(a >> 63)) & b);
         a <<= 1;
    }
    return c;
}

fq_elem fq_sqr(fq_elem a)
{
    return fq_mul(a, a);
}
#else
fq_elem fq_mul(fq_elem a, fq_elem b)
{
    fq_elem c;
    asm
    (
        "movq $0x807, %%r8\n"
        "pclmullqlqdq %2, %1\n"
        "movq %%r8, %%xmm1\n"
        "pclmullqhqdq %1, %%xmm1\n"
        "movq %%r8, %%xmm2\n"
        "pclmullqhqdq %%xmm1, %%xmm2\n"
        "movq %1, %0\n"
        "movq %%xmm1, %%r8\n"
        "xorq %%r8, %0\n"
        "movq %%xmm2, %%r8\n"
        "xorq %%r8, %0\n"
        : "=r" (c), "+x" (a)
        : "x"  (b)
        : "%r8", "%xmm1", "%xmm2"
    );
    return c;
}

fq_elem fq_sqr(fq_elem a)
{
    fq_elem c;
    asm
    (
        "movq $0x807, %%r8\n"
        "pclmullqlqdq %1, %1\n"
        "movq %%r8, %%xmm1\n"
        "pclmullqhqdq %1, %%xmm1\n"
        "movq %%r8, %%xmm2\n"
        "pclmullqhqdq %%xmm1, %%xmm2\n"
        "movq %1, %0\n"
        "movq %%xmm1, %%r8\n"
        "xorq %%r8, %0\n"
        "movq %%xmm2, %%r8\n"
        "xorq %%r8, %0\n"
        : "=r" (c), "+x" (a)
        :
        : "%r8", "%xmm1", "%xmm2"
    );
    return c;
}
#endif

fq_elem fq_mul_t(fq_elem a)
{
    return (a << 1) ^ ((-(a >> 63)) & fq_t64);
}

fq_elem fq_pow(fq_elem a, fq_expn n)
{
    fq_elem c;
    unsigned int i;
    c = fq_1;
    for (i=0; i<64; i++)
    {
        c = fq_sqr(c);
        if (n >> 63)
            c = fq_mul(a, c);
        n <<= 1;
    }
    return c;
}

fq_elem fq_inv(fq_elem a)
{
    return fq_pow(a, fq_qm2);
}

fq_elem fq_pow_2_k(fq_elem a, unsigned int k)
{
    while (k--)
        a = fq_mul(a, a);
    return a;
}

fq_elem fq_rnd(void)
{
    fq_elem a;
    randombytes((unsigned char*)&a, 8);
    return a;
}

void fq_print(fq_elem a)
{
    unsigned int i, first;
    printf("[");
    first = 1;
    for (i=0; i<64; i++)
    {
        if (a >> 63)
        {
            if (!first)
                printf("+");
            if (i < 62)
                printf("t^%d", 63-i);
            else if (i == 62)
                printf("t");
            else
                printf("1");
            first = 0;
        }
        a <<= 1;
    }
    if (first)
        printf("0");
    printf("]");
}

const unsigned char *fq_parse(fq_elem *x, const unsigned char *s)
{
    unsigned int i;
    if (!s)
        return NULL;
    *x = 0;
    for (i=0; i<sizeof(fq_elem); i++)
        *x += (fq_elem)(*s++) << (8*i);
    return s;
}

unsigned char *fq_serialize(unsigned char *s, fq_elem x)
{
    unsigned int i;
    for (i=0; i<sizeof(fq_elem); i++)
        *s++ = x >> (8*i);
    return s;
}

const fq2_elem fq2_0    = { fq_0, fq_0 };
const fq2_elem fq2_1    = { fq_1, fq_0 };
const fq2_elem fq2_u    = { fq_0, fq_1 };
const fq2_elem fq2_u2   = { fq_1, fq_t };  /* Fq2 = Fq[u]/<u^2+t*u+1> */
const fq2_expn fq2_q2m1 = { fq_qm1, fq_qm1 };
const fq2_expn fq2_q2m2 = { fq_qm2, fq_qm1 };
const fq2_expn fq2_inv_delta = {
    UINT64_C(0x0200080020008001),
    UINT64_C(0x0004001000400100)
};

void fq2_set(fq_elem *a, const fq_elem *b)
{
    a[0] = b[0];
    a[1] = b[1];
}

void fq2_add(fq_elem *a, const fq_elem *b, const fq_elem *c)
{
    a[0] = fq_add(b[0], c[0]);
    a[1] = fq_add(b[1], c[1]);
}

void fq2_mul(fq_elem *a, const fq_elem *b, const fq_elem *c)
{
    fq_elem d[3];
    d[0] = fq_mul(b[0], c[0]);
    d[1] = fq_add(fq_mul(b[0], c[1]), fq_mul(b[1], c[0]));
    d[2] = fq_mul(b[1], c[1]);
    a[0] = fq_add(d[0], /*fq_mul(fq2_u2[0], d[2])*/ d[2]); /* fq2_u2[0] = fq_1 */
    a[1] = fq_add(d[1], /*fq_mul(fq2_u2[1], d[2])*/ fq_mul_t(d[2])); /* fq2_u2[1] = fq_t */
}

void fq2_pow(fq_elem *a, const fq_elem *b, const fq_expn *c)
{
    fq2_elem d;
    fq_expn n;
    unsigned int i;
    fq2_set(d, fq2_1);
    n = c[1];
    for (i=0; i<64; i++)
    {
        fq2_mul(d, d, d);
        if (n >> 63)
            fq2_mul(d, d, b);
        n <<= 1;
    }
    n = c[0];
    for (i=0; i<64; i++)
    {
        fq2_mul(d, d, d);
        if (n >> 63)
            fq2_mul(d, d, b);
        n <<= 1;
    }
    fq2_set(a, d);
}

/*
void fq2_inv(fq_elem *a, const fq_elem *b)
{
    fq2_pow(a, b, fq2_q2m2);
}
*/

void fq2_inv(fq_elem *a, const fq_elem *b)
{
    fq_elem det1, det2, det3, det, inv;
    det1 = fq_add(b[0], fq_mul_t(b[1]) /*fq_mul(b[1], fq2_u2[1])*/);
    det2 = fq_mul(det1, b[0]);
    det3 = fq_mul(b[1], b[1]); /*fq_mul(fq_mul(b[1], b[1]), fq2_u2[0]);*/
    det  = fq_add(det2, det3);
    inv  = fq_inv(det);
    a[0] = fq_mul(det1, inv);
    a[1] = fq_mul(b[1], inv);
}

void fq2_pow_2_k(fq_elem *a, const fq_elem *b, unsigned int k)
{
    fq2_elem d;
    fq2_set(d, b);
    while (k--)
        fq2_mul(d, d, d);
    fq2_set(a, d);
}

void fq2_rnd(fq_elem *a)
{
    a[0] = fq_rnd();
    a[1] = fq_rnd();
}

void fq2_print(const fq_elem *a)
{
    printf("[");
    fq_print(a[1]);
    printf("*u+");
    fq_print(a[0]);
    printf("]");
}

const unsigned char *fq2_parse(fq_elem *x, const unsigned char *s)
{
    s = fq_parse(&x[0], s);
    s = fq_parse(&x[1], s);
    return s;
}

unsigned char *fq2_serialize(unsigned char *s, const fq_elem *x)
{
    s = fq_serialize(s, x[0]);
    s = fq_serialize(s, x[1]);
    return s;
}

int  fq_matrix_2x2_inv(fq_elem *inv, const fq_elem *mat)
{
    fq_elem dt, di, tmp;
    dt = fq_add(fq_mul(mat[0], mat[3]), fq_mul(mat[1], mat[2]));
    if (!dt) return 1;
    di = fq_inv(dt);
    tmp    = fq_mul(di, mat[3]);
    inv[1] = fq_mul(di, mat[1]);
    inv[2] = fq_mul(di, mat[2]);
    inv[3] = fq_mul(di, mat[0]);
    inv[0] = tmp;
    return 0;
}

void fq_matrix_2x2_map(fq_elem *x, const fq_elem *m, const fq_elem *y)
{
    fq_elem tmp1, tmp2;
    tmp1 = fq_add(fq_mul(m[0], y[0]), fq_mul(m[1], y[1]));
    tmp2 = fq_add(fq_mul(m[2], y[0]), fq_mul(m[3], y[1]));
    x[0] = tmp1;
    x[1] = tmp2;
}

void fq_matrix_2x2_rnd(fq_elem *mat, fq_elem *inv)
{
    do
    {
        mat[0] = fq_rnd();
        mat[1] = fq_rnd();
        mat[2] = fq_rnd();
        mat[3] = fq_rnd();
    }
    while (fq_matrix_2x2_inv(inv, mat));
}

const unsigned char *fq_matrix_2x2_parse(fq_elem *m, const unsigned char *s)
{
    s = fq_parse(&m[0], s);
    s = fq_parse(&m[1], s);
    s = fq_parse(&m[2], s);
    s = fq_parse(&m[3], s);
    return s;
}

unsigned char *fq_matrix_2x2_serialize(unsigned char *s, const fq_elem *m)
{
    s = fq_serialize(s, m[0]);
    s = fq_serialize(s, m[1]);
    s = fq_serialize(s, m[2]);
    s = fq_serialize(s, m[3]);
    return s;
}

void generate_skey(struct skey_t *sk)
{
    unsigned int i;
    memset(sk, 0, sizeof(struct skey_t));
    for (i=0; i<4; i++)
    {
        fq_matrix_2x2_rnd(&sk->L1[i][0][0], &sk->L1_inv[i][0][0]);
        fq_matrix_2x2_rnd(&sk->L2[i][0][0], &sk->L2_inv[i][0][0]);
        fq_matrix_2x2_rnd(&sk->L3[i][0][0], &sk->L3_inv[i][0][0]);
        fq_matrix_2x2_rnd(&sk->L4[i][0][0], &sk->L4_inv[i][0][0]);
    }
    fq2_rnd(sk->A3_0);
    randombytes(sk->a, 7);
    randombytes(sk->b, 8);
    randombytes(sk->c, 8);
    sk->a[0] &= 0x7f;
    sk->a[1] &= 0x7f;
    sk->a[3] &= 0x7f;
    sk->a[4] &= 0x7f;
    sk->a[5] &= 0x7f;
    sk->a[6] &= 0x7f;
    sk->b[0] &= 0x7f;
    sk->b[2] &= 0x7f;
    sk->b[4] &= 0x7f;
    sk->b[6] &= 0x7f;
    sk->c[0] &= 0x7f;
    sk->c[2] &= 0x7f;
    sk->c[4] &= 0x7f;
    sk->c[6] &= 0x7f;
    sk->a[2] = sk->a[0];
    sk->b[1] = sk->b[0];
    sk->b[3] = (sk->a[0] + sk->b[2] - sk->a[4]) & 0x7f;
    sk->b[5] = (sk->a[0] + sk->b[4] - sk->a[6]) & 0x7f;
    sk->b[7] = (sk->a[4] + sk->b[6] - sk->a[6]) & 0x7f;
    sk->c[1] = (sk->c[0] + sk->b[0] + sk->a[0] - sk->a[4] - sk->b[6] + delta) & 0x7f;
    sk->c[3] = (sk->b[0] + sk->c[2] - sk->b[4]) & 0x7f;
    sk->c[5] = (sk->b[2] + sk->c[4] - sk->b[4]) & 0x7f;
    sk->c[7] = (sk->a[0] + sk->b[2] + sk->c[6] - sk->a[4] - sk->b[6]) & 0x7f;
}

const unsigned char *skey_parse(struct skey_t *sk, const unsigned char *s)
{
    unsigned int i;
    memset(sk, 0, sizeof(struct skey_t));
    for (i=0; i<4; i++)
    {
        s = fq_matrix_2x2_parse(&sk->L1_inv[i][0][0], s);
        if (fq_matrix_2x2_inv(&sk->L1[i][0][0], &sk->L1_inv[i][0][0]))
            return NULL;
        s = fq_matrix_2x2_parse(&sk->L2_inv[i][0][0], s);
        if (fq_matrix_2x2_inv(&sk->L2[i][0][0], &sk->L2_inv[i][0][0]))
            return NULL;
        s = fq_matrix_2x2_parse(&sk->L3_inv[i][0][0], s);
        if (fq_matrix_2x2_inv(&sk->L3[i][0][0], &sk->L3_inv[i][0][0]))
            return NULL;
        s = fq_matrix_2x2_parse(&sk->L4_inv[i][0][0], s);
        if (fq_matrix_2x2_inv(&sk->L4[i][0][0], &sk->L4_inv[i][0][0]))
            return NULL;
    }
    s = fq2_parse(sk->A3_0, s);
    if ((sk->a[0] = *s++) >= 128) return NULL;
    if ((sk->a[1] = *s++) >= 128) return NULL;
    if ((sk->a[3] = *s++) >= 128) return NULL;
    if ((sk->a[4] = *s++) >= 128) return NULL;
    if ((sk->a[5] = *s++) >= 128) return NULL;
    if ((sk->a[6] = *s++) >= 128) return NULL;
    if ((sk->b[0] = *s++) >= 128) return NULL;
    if ((sk->b[2] = *s++) >= 128) return NULL;
    if ((sk->b[4] = *s++) >= 128) return NULL;
    if ((sk->b[6] = *s++) >= 128) return NULL;
    if ((sk->c[0] = *s++) >= 128) return NULL;
    if ((sk->c[2] = *s++) >= 128) return NULL;
    if ((sk->c[4] = *s++) >= 128) return NULL;
    if ((sk->c[6] = *s++) >= 128) return NULL;
    sk->a[2] = sk->a[0];
    sk->b[1] = sk->b[0];
    sk->b[3] = (sk->a[0] + sk->b[2] - sk->a[4]) & 0x7f;
    sk->b[5] = (sk->a[0] + sk->b[4] - sk->a[6]) & 0x7f;
    sk->b[7] = (sk->a[4] + sk->b[6] - sk->a[6]) & 0x7f;
    sk->c[1] = (sk->c[0] + sk->b[0] + sk->a[0] - sk->a[4] - sk->b[6] + delta) & 0x7f;
    sk->c[3] = (sk->b[0] + sk->c[2] - sk->b[4]) & 0x7f;
    sk->c[5] = (sk->b[2] + sk->c[4] - sk->b[4]) & 0x7f;
    sk->c[7] = (sk->a[0] + sk->b[2] + sk->c[6] - sk->a[4] - sk->b[6]) & 0x7f;
    return s;
}

unsigned char *skey_serialize(unsigned char *s, const struct skey_t *sk)
{
    unsigned int i;
    for (i=0; i<4; i++)
    {
        s = fq_matrix_2x2_serialize(s, &sk->L1_inv[i][0][0]);
        s = fq_matrix_2x2_serialize(s, &sk->L2_inv[i][0][0]);
        s = fq_matrix_2x2_serialize(s, &sk->L3_inv[i][0][0]);
        s = fq_matrix_2x2_serialize(s, &sk->L4_inv[i][0][0]);
    }
    s = fq2_serialize(s, sk->A3_0);
    *s++ = sk->a[0];
    *s++ = sk->a[1];
    *s++ = sk->a[3];
    *s++ = sk->a[4];
    *s++ = sk->a[5];
    *s++ = sk->a[6];
    *s++ = sk->b[0];
    *s++ = sk->b[2];
    *s++ = sk->b[4];
    *s++ = sk->b[6];
    *s++ = sk->c[0];
    *s++ = sk->c[2];
    *s++ = sk->c[4];
    *s++ = sk->c[6];
    return s; 
}

void apply_expn_1(fq_elem *x, const unsigned char *a, const fq_elem *y)
{
    fq2_elem tmp1, tmp2;
    fq2_pow_2_k(tmp1, &y[0], a[0]);
    fq2_pow_2_k(tmp2, &y[2], a[1]);
    fq2_mul(&x[0], tmp1, tmp2);
    fq2_pow_2_k(tmp1, &y[0], a[2]);
    fq2_pow_2_k(tmp2, &y[4], a[3]);
    fq2_mul(&x[2], tmp1, tmp2);
    fq2_pow_2_k(tmp1, &y[0], a[4]);
    fq2_pow_2_k(tmp2, &y[6], a[5]);
    fq2_mul(&x[4], tmp1, tmp2);
    fq2_pow_2_k(&x[6], &y[0], a[6]);
}

void apply_expn_2(fq_elem *x, const unsigned char *b, const fq_elem *y)
{
    fq2_elem tmp1, tmp2;
    fq2_pow_2_k(tmp1, &y[0], b[0]);
    fq2_pow_2_k(tmp2, &y[2], b[1]);
    fq2_mul(&x[0], tmp1, tmp2);
    fq2_pow_2_k(tmp1, &y[0], b[2]);
    fq2_pow_2_k(tmp2, &y[4], b[3]);
    fq2_mul(&x[2], tmp1, tmp2);
    fq2_pow_2_k(tmp1, &y[0], b[4]);
    fq2_pow_2_k(tmp2, &y[6], b[5]);
    fq2_mul(&x[4], tmp1, tmp2);
    fq2_pow_2_k(tmp1, &y[4], b[6]);
    fq2_pow_2_k(tmp2, &y[6], b[7]);
    fq2_mul(&x[6], tmp1, tmp2);
}

void apply_expn_3(fq_elem *x, const unsigned char *c, const fq_elem *y)
{
    fq2_elem tmp1, tmp2;
    fq2_pow_2_k(tmp1, &y[0], c[0]);
    fq2_pow_2_k(tmp2, &y[6], c[1]);
    fq2_mul(&x[0], tmp1, tmp2);
    fq2_pow_2_k(tmp1, &y[0], c[2]);
    fq2_pow_2_k(tmp2, &y[4], c[3]);
    fq2_mul(&x[2], tmp1, tmp2);
    fq2_pow_2_k(tmp1, &y[2], c[4]);
    fq2_pow_2_k(tmp2, &y[4], c[5]);
    fq2_mul(&x[4], tmp1, tmp2);
    fq2_pow_2_k(tmp1, &y[2], c[6]);
    fq2_pow_2_k(tmp2, &y[6], c[7]);
    fq2_mul(&x[6], tmp1, tmp2);
}

void encrypt_with_skey(fq_elem *ct, const struct skey_t *sk, const fq_elem *pt)
{
    unsigned int i;
    fq_elem t1[8], t2[8];
    for (i=0; i<4; i++)
        fq_matrix_2x2_map(&t1[2*i], &sk->L1[i][0][0], &pt[2*i]);
    apply_expn_1(t2, sk->a, t1);
    for (i=0; i<4; i++)
        fq_matrix_2x2_map(&t1[2*i], &sk->L2[i][0][0], &t2[2*i]);
    apply_expn_2(t2, sk->b, t1);
    for (i=0; i<4; i++)
        fq_matrix_2x2_map(&t1[2*i], &sk->L3[i][0][0], &t2[2*i]);
    fq2_add(&t1[0], &t1[0], sk->A3_0);
    apply_expn_3(t2, sk->c, t1);
    for (i=0; i<4; i++)
        fq_matrix_2x2_map(&ct[2*i], &sk->L4[i][0][0], &t2[2*i]);
}

void apply_expn_inv_1(fq_elem *x, const unsigned char *a, const fq_elem *y)
{
    fq2_elem tmp1, tmp2, inv1;
    fq2_inv(inv1, &y[6]);
    fq2_pow_2_k(&x[0], &y[6], -a[6] & 0x7f);
    fq2_pow_2_k(tmp1,  &y[0], -a[1] & 0x7f);
    fq2_pow_2_k(tmp2,  inv1, (a[0]-a[1]-a[6]) & 0x7f);
    fq2_mul(&x[2], tmp1, tmp2);
    fq2_pow_2_k(tmp1,  &y[2], -a[3] & 0x7f);
    fq2_pow_2_k(tmp2,  inv1, (a[2]-a[3]-a[6]) & 0x7f);
    fq2_mul(&x[4], tmp1, tmp2);
    fq2_pow_2_k(tmp1,  &y[4], -a[5] & 0x7f);
    fq2_pow_2_k(tmp2,  inv1, (a[4]-a[5]-a[6]) & 0x7f);
    fq2_mul(&x[6], tmp1, tmp2);
}

void apply_expn_inv_2(fq_elem *x, const unsigned char *b, const fq_elem *y)
{
    fq2_elem tmp1, tmp2, tmp3, tmp4, inv1, inv2, inv3;
    fq2_inv(inv1, &y[2]);
    fq2_inv(inv2, &y[4]);
    fq2_inv(inv3, &y[6]);
    fq2_pow_2_k(tmp1, &y[2], (-1-b[2]) & 0x7f);
    fq2_pow_2_k(tmp2, &y[4], (-1-b[4]) & 0x7f);
    fq2_pow_2_k(tmp3, inv3,  (-1-b[4]+b[5]-b[7]) & 0x7f);
    fq2_mul(tmp1, tmp1, tmp2);
    fq2_mul(&x[0], tmp1, tmp3);
    fq2_pow_2_k(tmp1, &y[0], -b[0] & 0x7f);
    fq2_pow_2_k(tmp2, inv1,  (-1-b[2]) & 0x7f);
    fq2_pow_2_k(tmp3, inv2,  (-1-b[4]) & 0x7f);
    fq2_pow_2_k(tmp4, &y[6], (-1-b[4]+b[5]-b[7]) & 0x7f);
    fq2_mul(tmp1, tmp1, tmp2);
    fq2_mul(tmp1, tmp1, tmp3);
    fq2_mul(&x[2], tmp1, tmp4);
    fq2_pow_2_k(tmp1, &y[2], (-1-b[3]) & 0x7f);
    fq2_pow_2_k(tmp2, inv2,  (-1-b[5]-b[6]+b[7]) & 0x7f);
    fq2_pow_2_k(tmp3, &y[6], (-1-b[6]) & 0x7f);
    fq2_mul(tmp1, tmp1, tmp2);
    fq2_mul(&x[4], tmp1, tmp3);
    fq2_pow_2_k(tmp1, inv1,  (-1-b[2]+b[4]-b[5]) & 0x7f);
    fq2_pow_2_k(tmp2, &y[4], (-1-b[5]) & 0x7f);
    fq2_pow_2_k(tmp3, &y[6], (-1-b[7]) & 0x7f);
    fq2_mul(tmp1, tmp1, tmp2);
    fq2_mul(&x[6], tmp1, tmp3);
}

void apply_expn_inv_3(fq_elem *x, const unsigned char *c, const fq_elem *y)
{
    fq2_elem tmp1, tmp2, tmp3, tmp4, inv1, inv2, inv3, inv4;
    fq2_inv(inv1, &y[0]);
    fq2_inv(inv2, &y[2]);
    fq2_inv(inv3, &y[4]);
    fq2_inv(inv4, &y[6]);
    fq2_pow_2_k(tmp1, inv1,  -c[0] & 0x7f);
    fq2_pow_2_k(tmp2, &y[2], (delta-c[2]) & 0x7f);
    fq2_pow_2_k(tmp3, inv3,  (delta-c[2]+c[3]-c[5]) & 0x7f);
    fq2_pow_2_k(tmp4, &y[6], (c[1]-c[7]-c[0]) & 0x7f);
    fq2_mul(tmp1, tmp1, tmp2);
    fq2_mul(tmp1, tmp1, tmp3);
    fq2_mul(tmp1, tmp1, tmp4);
    fq2_pow(&x[0], tmp1, fq2_inv_delta);
    fq2_pow_2_k(tmp1, inv1,  (delta-c[1]+c[7]-c[6]) & 0x7f);
    fq2_pow_2_k(tmp2, &y[2], (delta-c[1]+c[7]+c[0]-c[6]-c[2]) & 0x7f);
    fq2_pow_2_k(tmp3, inv3,  -c[4] & 0x7f);
    fq2_pow_2_k(tmp4, &y[6], (delta-c[6]) & 0x7f);
    fq2_mul(tmp1, tmp1, tmp2);
    fq2_mul(tmp1, tmp1, tmp3);
    fq2_mul(tmp1, tmp1, tmp4);
    fq2_pow(&x[2], tmp1, fq2_inv_delta);
    fq2_pow_2_k(tmp1, &y[0], (delta-c[1]+c[7]-c[5]+c[4]-c[6]) & 0x7f);
    fq2_pow_2_k(tmp2, inv2,  (delta-c[1]+c[7]-c[5]+c[4]-c[6]+c[0]-c[2]) & 0x7f);
    fq2_pow_2_k(tmp3, &y[4], (delta-c[5]) & 0x7f);
    fq2_pow_2_k(tmp4, inv4,  (delta-c[5]+c[4]-c[6]) & 0x7f);
    fq2_mul(tmp1, tmp1, tmp2);
    fq2_mul(tmp1, tmp1, tmp3);
    fq2_mul(tmp1, tmp1, tmp4);
    fq2_pow(&x[4], tmp1, fq2_inv_delta);
    fq2_pow_2_k(tmp1, &y[0], (delta-c[1]) & 0x7f);
    fq2_pow_2_k(tmp2, inv2,  (delta-c[1]+c[0]-c[2]) & 0x7f);
    fq2_pow_2_k(tmp3, &y[4], (delta-c[1]-c[5]+c[3]+c[0]-c[2]) & 0x7f);
    fq2_pow_2_k(tmp4, inv4,  (delta-c[1]-c[5]+c[3]+c[0]-c[2]+c[4]-c[6]) & 0x7f);
    fq2_mul(tmp1, tmp1, tmp2);
    fq2_mul(tmp1, tmp1, tmp3);
    fq2_mul(tmp1, tmp1, tmp4);
    fq2_pow(&x[6], tmp1, fq2_inv_delta);
}

void decrypt_with_skey(fq_elem *pt, const struct skey_t *sk, const fq_elem *ct)
{
    unsigned int i;
    fq_elem t1[8], t2[8];
    for (i=0; i<4; i++)
        fq_matrix_2x2_map(&t2[2*i], &sk->L4_inv[i][0][0], &ct[2*i]);
    apply_expn_inv_3(t1, sk->c, t2);
    fq2_add(&t1[0], &t1[0], sk->A3_0);
    for (i=0; i<4; i++)
        fq_matrix_2x2_map(&t2[2*i], &sk->L3_inv[i][0][0], &t1[2*i]);
    apply_expn_inv_2(t1, sk->b, t2);
    for (i=0; i<4; i++)
        fq_matrix_2x2_map(&t2[2*i], &sk->L2_inv[i][0][0], &t1[2*i]);
    apply_expn_inv_1(t1, sk->a, t2);
    for (i=0; i<4; i++)
        fq_matrix_2x2_map(&pt[2*i], &sk->L1_inv[i][0][0], &t1[2*i]);
}

void fq2_matrix_2x2_map_polynomial(fq_elem *p, const fq_elem *m, const fq_elem *q,
    unsigned int n)
{
    unsigned int i;
    for (i=0; i<n; i++)
        fq_matrix_2x2_map(&p[2*i], m, &q[2*i]);
}

void fq2_polynomial_product(fq_elem *p,
    const fq_elem *q1, unsigned int a1, unsigned int n1,
    const fq_elem *q2, unsigned int a2, unsigned int n2)
{
    unsigned int i, j, k;
    fq2_elem tmp1, tmp2;
    for (i=k=0; i<n1; i++)
        for (j=0; j<n2; j++, k++)
        {
            fq2_pow_2_k(tmp1, &q1[2*i], a1);
            fq2_pow_2_k(tmp2, &q2[2*j], a2);
            fq2_mul(&p[2*k], tmp1, tmp2);
        }
}

void fq2_polynomial_pow_2_k(fq_elem *p,
    const fq_elem *q, unsigned int k, unsigned int n)
{
    unsigned int i;
    for (i=0; i<n; i++)
        fq2_pow_2_k(&p[2*i], &q[2*i], k);
}

const unsigned int p12_3_red_idx[16] =
{
    0, 1, 2, 3, 4, 5, 6, 7, 2, 3, 8, 9, 6, 7, 10, 11
};

const unsigned int p34_3_red_idx[16] =
{
    0, 1, 2, 3, 4, 5, 6, 7, 2, 3, 8, 9, 6, 7, 10, 11
};

const unsigned int p56_3_red_idx[8] =
{
    0, 1, 2, 3, 1, 4, 3, 5
};

const unsigned int p78_3_red_idx[8] =
{
    0, 1, 2, 3, 1, 4, 3, 5
};

const unsigned int p12_4_red_idx[78] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
    65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77
};

const unsigned int p34_4_red_idx[78] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,  1,
     4,  3,  5, 12, 13,  7, 10,  9, 11, 14, 15,  2,  3,
    16, 17,  5, 18,  8,  9, 19, 20, 11, 21,  3,  5, 17,
    18, 13, 22,  9, 11, 20, 21, 15, 23,  4, 12,  5, 13,
    24, 25, 10, 14, 11, 15, 26, 27,  5, 13, 18, 22, 25,
    28, 11, 15, 21, 23, 27, 29, 30, 31, 32, 33, 34, 35
};

const unsigned int p56_4_red_idx[72] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
     1,  4,  3,  5, 12, 13,  7, 10,  9, 11, 14, 15,
     2,  3, 16, 17,  5, 18,  8,  9, 19, 20, 11, 21,
     3,  5, 17, 18, 13, 22,  9, 11, 20, 21, 15, 23,
     4, 12,  5, 13, 24, 25, 10, 14, 11, 15, 26, 27,
     5, 13, 18, 22, 25, 28, 11, 15, 21, 23, 27, 29
};

const unsigned int p78_4_red_idx[72] =
{
     0,  1,  2,  3,  4,  5,  2,  3,  6,  7,  5,  8,
     1,  4,  3,  5,  9, 10,  3,  5,  7,  8, 10, 11,
    12, 13, 14, 15, 16, 17, 14, 15, 18, 19, 17, 20,
    13, 16, 15, 17, 21, 22, 15, 17, 19, 20, 22, 23,
     4,  9,  5, 10, 24, 25,  5, 10,  8, 11, 25, 26,
    16, 21, 17, 22, 27, 28, 17, 22, 20, 23, 28, 29
};

void generate_pkey(struct pkey_t *pk, struct skey_t *sk)
{
    unsigned int i;
    fq2_elem p12_1[2],  p34_1[2], p56_1[2], p78_1[2];
    fq2_elem p12_2[4],  p34_2[4], p56_2[4], p78_2[2];
    fq2_elem p12_3[16], p34_3[16], p56_3[8], p78_3[8];
    fq2_elem p12_3_red[13], p34_3_red[12], p56_3_red[6], p78_3_red[6];
    fq2_elem p12_4[78], p34_4[78], p56_4[72], p78_4[72];
    memset(pk, 0, sizeof(struct pkey_t));
    fq2_set(p12_1[0], fq2_1);
    fq2_set(p12_1[1], fq2_u);
    fq2_set(p34_1[0], fq2_1);
    fq2_set(p34_1[1], fq2_u);
    fq2_set(p56_1[0], fq2_1);
    fq2_set(p56_1[1], fq2_u);
    fq2_set(p78_1[0], fq2_1);
    fq2_set(p78_1[1], fq2_u);
    fq2_matrix_2x2_map_polynomial(p12_1[0], &sk->L1[0][0][0], p12_1[0], 2);
    fq2_matrix_2x2_map_polynomial(p34_1[0], &sk->L1[1][0][0], p34_1[0], 2);
    fq2_matrix_2x2_map_polynomial(p56_1[0], &sk->L1[2][0][0], p56_1[0], 2);
    fq2_matrix_2x2_map_polynomial(p78_1[0], &sk->L1[3][0][0], p78_1[0], 2);
    fq2_polynomial_product(p12_2[0], p12_1[0], sk->a[0], 2, p34_1[0], sk->a[1], 2);
    fq2_polynomial_product(p34_2[0], p12_1[0], sk->a[2], 2, p56_1[0], sk->a[3], 2);
    fq2_polynomial_product(p56_2[0], p12_1[0], sk->a[4], 2, p78_1[0], sk->a[5], 2);
    fq2_polynomial_pow_2_k(p78_2[0], p12_1[0], sk->a[6], 2);
    fq2_matrix_2x2_map_polynomial(p12_2[0], &sk->L2[0][0][0], p12_2[0], 4);
    fq2_matrix_2x2_map_polynomial(p34_2[0], &sk->L2[1][0][0], p34_2[0], 4);
    fq2_matrix_2x2_map_polynomial(p56_2[0], &sk->L2[2][0][0], p56_2[0], 4);
    fq2_matrix_2x2_map_polynomial(p78_2[0], &sk->L2[3][0][0], p78_2[0], 2);
    fq2_polynomial_product(p12_3[0], p12_2[0], sk->b[0], 4, p34_2[0], sk->b[1], 4);
    fq2_polynomial_product(p34_3[0], p12_2[0], sk->b[2], 4, p56_2[0], sk->b[3], 4);
    fq2_polynomial_product(p56_3[0], p12_2[0], sk->b[4], 4, p78_2[0], sk->b[5], 2);
    fq2_polynomial_product(p78_3[0], p56_2[0], sk->b[6], 4, p78_2[0], sk->b[7], 2);
    memset(p12_3_red, 0, sizeof(p12_3_red));
    memset(p34_3_red, 0, sizeof(p34_3_red));
    memset(p56_3_red, 0, sizeof(p56_3_red));
    memset(p78_3_red, 0, sizeof(p78_3_red));
    for (i=0; i<16; i++)
        fq2_add(p12_3_red[p12_3_red_idx[i]], p12_3_red[p12_3_red_idx[i]], p12_3[i]);
    for (i=0; i<16; i++)
        fq2_add(p34_3_red[p34_3_red_idx[i]], p34_3_red[p34_3_red_idx[i]], p34_3[i]);
    for (i=0; i<8; i++)
        fq2_add(p56_3_red[p56_3_red_idx[i]], p56_3_red[p56_3_red_idx[i]], p56_3[i]);
    for (i=0; i<8; i++)
        fq2_add(p78_3_red[p78_3_red_idx[i]], p78_3_red[p78_3_red_idx[i]], p78_3[i]);
    fq2_matrix_2x2_map_polynomial(p12_3_red[0], &sk->L3[0][0][0], p12_3_red[0], 12);
    fq2_matrix_2x2_map_polynomial(p34_3_red[0], &sk->L3[1][0][0], p34_3_red[0], 12);
    fq2_matrix_2x2_map_polynomial(p56_3_red[0], &sk->L3[2][0][0], p56_3_red[0], 6);
    fq2_matrix_2x2_map_polynomial(p78_3_red[0], &sk->L3[3][0][0], p78_3_red[0], 6);
    fq2_set(p12_3_red[12], sk->A3_0);
    fq2_polynomial_product(p12_4[0], p12_3_red[0], sk->c[0], 13, p78_3_red[0], sk->c[1], 6);
    fq2_polynomial_product(p34_4[0], p12_3_red[0], sk->c[2], 13, p56_3_red[0], sk->c[3], 6);
    fq2_polynomial_product(p56_4[0], p34_3_red[0], sk->c[4], 12, p56_3_red[0], sk->c[5], 6);
    fq2_polynomial_product(p78_4[0], p34_3_red[0], sk->c[6], 12, p78_3_red[0], sk->c[7], 6);
    for (i=0; i<78; i++)
        fq2_add(pk->p12[p12_4_red_idx[i]], pk->p12[p12_4_red_idx[i]], p12_4[i]);
    for (i=0; i<78; i++)
        fq2_add(pk->p34[p34_4_red_idx[i]], pk->p34[p34_4_red_idx[i]], p34_4[i]);
    for (i=0; i<72; i++)
        fq2_add(pk->p56[p56_4_red_idx[i]], pk->p56[p56_4_red_idx[i]], p56_4[i]);
    for (i=0; i<72; i++)
        fq2_add(pk->p78[p78_4_red_idx[i]], pk->p78[p78_4_red_idx[i]], p78_4[i]);
    fq2_matrix_2x2_map_polynomial(pk->p12[0], &sk->L4[0][0][0], pk->p12[0], 78);
    fq2_matrix_2x2_map_polynomial(pk->p34[0], &sk->L4[1][0][0], pk->p34[0], 36);
    fq2_matrix_2x2_map_polynomial(pk->p56[0], &sk->L4[2][0][0], pk->p56[0], 30);
    fq2_matrix_2x2_map_polynomial(pk->p78[0], &sk->L4[3][0][0], pk->p78[0], 30);
    /* x1, x2 */
    pk->f[0]  = (sk->a[4] + sk->b[6] + sk->c[1]) & 0x3f;
    pk->f[1]  = (sk->a[0] + sk->b[0] + sk->c[0]) & 0x3f;
    pk->f[2]  = (sk->a[0] + sk->b[0] + sk->c[2]) & 0x3f;
    pk->f[3]  = (sk->a[0] + sk->b[2] + sk->c[4]) & 0x3f;
    pk->f[4]  = (sk->a[0] + sk->b[2] + sk->c[6]) & 0x3f;
    /* x3, x4 */
    pk->f[5]  = (sk->a[1] + sk->b[0] + sk->c[0]) & 0x3f;
    pk->f[6]  = (sk->a[1] + sk->b[0] + sk->c[2]) & 0x3f;
    pk->f[7]  = (sk->a[1] + sk->b[2] + sk->c[4]) & 0x3f;
    pk->f[8]  = (sk->a[1] + sk->b[2] + sk->c[6]) & 0x3f;
    /* x5, x6 */
    pk->f[9]  = (sk->a[3] + sk->b[0] + sk->c[0]) & 0x3f;
    pk->f[10] = (sk->a[3] + sk->b[0] + sk->c[2]) & 0x3f;
    /* x7, x8 */
    pk->f[11] = (sk->a[5] + sk->b[6] + sk->c[1]) & 0x3f;
    pk->f[12] = (sk->a[5] + sk->b[3] + sk->c[4]) & 0x3f;
    pk->f[13] = (sk->a[5] + sk->b[3] + sk->c[6]) & 0x3f;
}

fq_elem fq_monomio(unsigned int n, ...)
{
    unsigned int i;
    fq_elem x, y;
    va_list ap;
    if (!n) return fq_1;
    va_start(ap, n);
    x = va_arg(ap, fq_elem);
    for (i=1; i<n; i++)
    {
        y = va_arg(ap, fq_elem);
        x = fq_mul(x, y);
    }
    va_end(ap);
    return x;
}

void encrypt_with_pkey(fq_elem *ct, const struct pkey_t *pk, const fq_elem *pt)
{
    unsigned int i;
    fq_elem y1[28], y2[28], y3[28], y4[28], c12[78], c34[36], c56[30], c78[30];
    for (i=0; i<5; i++)
    {
        y1[i]    = fq_pow_2_k(pt[0], pk->f[i]);
        y1[i+5]  = fq_pow_2_k(pt[1], pk->f[i]);
    }
    for (i=0; i<4; i++)
    {
        y1[i+10] = fq_pow_2_k(pt[2], pk->f[i+5]);
        y1[i+14] = fq_pow_2_k(pt[3], pk->f[i+5]);
    }
    for (i=0; i<2; i++)
    {
        y1[i+18] = fq_pow_2_k(pt[4], pk->f[i+9]);
        y1[i+20] = fq_pow_2_k(pt[5], pk->f[i+9]);
    }
    for (i=0; i<3; i++)
    {
        y1[i+22] = fq_pow_2_k(pt[6], pk->f[i+11]);
        y1[i+25] = fq_pow_2_k(pt[7], pk->f[i+11]);
    }    
    for (i=0; i<28; i++)
    {
        y2[i] = fq_mul(y1[i], y1[i]);
        y3[i] = fq_mul(y1[i], y2[i]);
        y4[i] = fq_mul(y1[i], y3[i]);
    }
    c12[0]  = fq_monomio(5, y2[0], y2[1], y1[22], y1[18], y1[10]);
    c12[1]  = fq_monomio(6, y1[5], y1[22], y1[0], y2[1], y1[18], y1[10]);
    c12[2]  = fq_monomio(5, y2[0], y2[1], y1[25], y1[18], y1[10]);
    c12[3]  = fq_monomio(6, y1[5], y1[25], y1[0], y2[1], y1[18], y1[10]);
    c12[4]  = fq_monomio(5, y2[5], y1[22], y1[18], y2[1], y1[10]);
    c12[5]  = fq_monomio(5, y2[5], y1[25], y1[18], y2[1], y1[10]);
    c12[6]  = fq_monomio(5, y2[0], y2[1], y1[22], y1[20], y1[10]);
    c12[7]  = fq_monomio(6, y1[5], y1[22], y1[0], y2[1], y1[20], y1[10]);
    c12[8]  = fq_monomio(5, y2[0], y2[1], y1[25], y1[20], y1[10]);
    c12[9]  = fq_monomio(6, y1[5], y1[25], y1[0], y2[1], y1[20], y1[10]);
    c12[10] = fq_monomio(5, y2[5], y1[22], y1[20], y2[1], y1[10]);
    c12[11] = fq_monomio(5, y2[5], y1[25], y1[20], y2[1], y1[10]);
    c12[12] = fq_monomio(6, y2[0], y1[1], y1[22], y1[18], y1[6], y1[10]);
    c12[13] = fq_monomio(7, y1[5], y1[6], y1[22], y1[0], y1[1], y1[18], y1[10]);
    c12[14] = fq_monomio(6, y2[0], y1[1], y1[25], y1[18], y1[6], y1[10]);
    c12[15] = fq_monomio(7, y1[5], y1[6], y1[25], y1[0], y1[1], y1[18], y1[10]);
    c12[16] = fq_monomio(6, y2[5], y1[6], y1[22], y1[18], y1[10], y1[1]);
    c12[17] = fq_monomio(6, y2[5], y1[6], y1[25], y1[18], y1[10], y1[1]);
    c12[18] = fq_monomio(6, y2[0], y1[1], y1[22], y1[20], y1[6], y1[10]);
    c12[19] = fq_monomio(7, y1[5], y1[6], y1[22], y1[0], y1[1], y1[20], y1[10]);
    c12[20] = fq_monomio(6, y2[0], y1[1], y1[25], y1[20], y1[6], y1[10]);
    c12[21] = fq_monomio(7, y1[5], y1[6], y1[25], y1[0], y1[1], y1[20], y1[10]);
    c12[22] = fq_monomio(6, y2[5], y1[6], y1[22], y1[20], y1[10], y1[1]);
    c12[23] = fq_monomio(6, y2[5], y1[6], y1[25], y1[20], y1[10], y1[1]);
    c12[24] = fq_monomio(5, y2[0], y2[1], y1[22], y1[18], y1[14]);
    c12[25] = fq_monomio(6, y1[5], y1[22], y1[0], y2[1], y1[18], y1[14]);
    c12[26] = fq_monomio(5, y2[0], y2[1], y1[25], y1[18], y1[14]);
    c12[27] = fq_monomio(6, y1[5], y1[25], y1[0], y2[1], y1[18], y1[14]);
    c12[28] = fq_monomio(5, y2[5], y1[22], y1[18], y2[1], y1[14]);
    c12[29] = fq_monomio(5, y2[5], y1[25], y1[18], y2[1], y1[14]);
    c12[30] = fq_monomio(5, y2[0], y2[1], y1[22], y1[20], y1[14]);
    c12[31] = fq_monomio(6, y1[5], y1[22], y1[0], y2[1], y1[20], y1[14]);
    c12[32] = fq_monomio(5, y2[0], y2[1], y1[25], y1[20], y1[14]);
    c12[33] = fq_monomio(6, y1[5], y1[25], y1[0], y2[1], y1[20], y1[14]);
    c12[34] = fq_monomio(5, y2[5], y1[22], y1[20], y2[1], y1[14]);
    c12[35] = fq_monomio(5, y2[5], y1[25], y1[20], y2[1], y1[14]);
    c12[36] = fq_monomio(6, y2[0], y1[1], y1[22], y1[18], y1[6], y1[14]);
    c12[37] = fq_monomio(7, y1[5], y1[6], y1[22], y1[0], y1[1], y1[18], y1[14]);
    c12[38] = fq_monomio(6, y2[0], y1[1], y1[25], y1[18], y1[6], y1[14]);
    c12[39] = fq_monomio(7, y1[5], y1[6], y1[25], y1[0], y1[1], y1[18], y1[14]);
    c12[40] = fq_monomio(6, y2[5], y1[6], y1[22], y1[18], y1[14], y1[1]);
    c12[41] = fq_monomio(6, y2[5], y1[6], y1[25], y1[18], y1[14], y1[1]);
    c12[42] = fq_monomio(6, y2[0], y1[1], y1[22], y1[20], y1[6], y1[14]);
    c12[43] = fq_monomio(7, y1[5], y1[6], y1[22], y1[0], y1[1], y1[20], y1[14]);
    c12[44] = fq_monomio(6, y2[0], y1[1], y1[25], y1[20], y1[6], y1[14]);
    c12[45] = fq_monomio(7, y1[5], y1[6], y1[25], y1[0], y1[1], y1[20], y1[14]);
    c12[46] = fq_monomio(6, y2[5], y1[6], y1[22], y1[20], y1[14], y1[1]);
    c12[47] = fq_monomio(6, y2[5], y1[6], y1[25], y1[20], y1[14], y1[1]);
    c12[48] = fq_monomio(5, y2[0], y1[22], y1[18], y2[6], y1[10]);
    c12[49] = fq_monomio(6, y1[5], y2[6], y1[22], y1[0], y1[18], y1[10]);
    c12[50] = fq_monomio(5, y2[0], y1[25], y1[18], y2[6], y1[10]);
    c12[51] = fq_monomio(6, y1[5], y2[6], y1[25], y1[0], y1[18], y1[10]);
    c12[52] = fq_monomio(5, y2[5], y2[6], y1[22], y1[18], y1[10]);
    c12[53] = fq_monomio(5, y2[5], y2[6], y1[25], y1[18], y1[10]);
    c12[54] = fq_monomio(5, y2[0], y1[22], y1[20], y2[6], y1[10]);
    c12[55] = fq_monomio(6, y1[5], y2[6], y1[22], y1[0], y1[20], y1[10]);
    c12[56] = fq_monomio(5, y2[0], y1[25], y1[20], y2[6], y1[10]);
    c12[57] = fq_monomio(6, y1[5], y2[6], y1[25], y1[0], y1[20], y1[10]);
    c12[58] = fq_monomio(5, y2[5], y2[6], y1[22], y1[20], y1[10]);
    c12[59] = fq_monomio(5, y2[5], y2[6], y1[25], y1[20], y1[10]);
    c12[60] = fq_monomio(5, y2[0], y1[22], y1[18], y2[6], y1[14]);
    c12[61] = fq_monomio(6, y1[5], y2[6], y1[22], y1[0], y1[18], y1[14]);
    c12[62] = fq_monomio(5, y2[0], y1[25], y1[18], y2[6], y1[14]);
    c12[63] = fq_monomio(6, y1[5], y2[6], y1[25], y1[0], y1[18], y1[14]);
    c12[64] = fq_monomio(5, y2[5], y2[6], y1[22], y1[18], y1[14]);
    c12[65] = fq_monomio(5, y2[5], y2[6], y1[25], y1[18], y1[14]);
    c12[66] = fq_monomio(5, y2[0], y1[22], y1[20], y2[6], y1[14]);
    c12[67] = fq_monomio(6, y1[5], y2[6], y1[22], y1[0], y1[20], y1[14]);
    c12[68] = fq_monomio(5, y2[0], y1[25], y1[20], y2[6], y1[14]);
    c12[69] = fq_monomio(6, y1[5], y2[6], y1[25], y1[0], y1[20], y1[14]);
    c12[70] = fq_monomio(5, y2[5], y2[6], y1[22], y1[20], y1[14]);
    c12[71] = fq_monomio(5, y2[5], y2[6], y1[25], y1[20], y1[14]);
    c12[72] = fq_monomio(2, y2[0], y1[22]);
    c12[73] = fq_monomio(3, y1[5], y1[22], y1[0]);
    c12[74] = fq_monomio(2, y2[0], y1[25]);
    c12[75] = fq_monomio(3, y1[5], y1[25], y1[0]);
    c12[76] = fq_monomio(2, y2[5], y1[22]);
    c12[77] = fq_monomio(2, y2[5], y1[25]);
    c34[0]  = fq_monomio(3, y4[2], y2[11], y1[19]);
    c34[1]  = fq_monomio(4, y1[7], y2[11], y3[2], y1[19]);
    c34[2]  = fq_monomio(4, y4[2], y1[15], y1[19], y1[11]);
    c34[3]  = fq_monomio(5, y1[7], y1[15], y3[2], y1[19], y1[11]);
    c34[4]  = fq_monomio(4, y2[7], y2[11], y1[19], y2[2]);
    c34[5]  = fq_monomio(5, y2[7], y1[15], y1[19], y2[2], y1[11]);
    c34[6]  = fq_monomio(3, y4[2], y2[11], y1[21]);
    c34[7]  = fq_monomio(4, y1[7], y2[11], y3[2], y1[21]);
    c34[8]  = fq_monomio(4, y4[2], y1[15], y1[21], y1[11]);
    c34[9]  = fq_monomio(5, y1[7], y1[15], y3[2], y1[21], y1[11]);
    c34[10] = fq_monomio(4, y2[7], y2[11], y1[21], y2[2]);
    c34[11] = fq_monomio(5, y2[7], y1[15], y1[21], y2[2], y1[11]);
    c34[12] = fq_monomio(4, y3[7], y2[11], y1[19], y1[2]);
    c34[13] = fq_monomio(5, y3[7], y1[15], y1[19], y1[11], y1[2]);
    c34[14] = fq_monomio(4, y3[7], y2[11], y1[21], y1[2]);
    c34[15] = fq_monomio(5, y3[7], y1[15], y1[21], y1[11], y1[2]);
    c34[16] = fq_monomio(3, y4[2], y2[15], y1[19]);
    c34[17] = fq_monomio(4, y1[7], y2[15], y3[2], y1[19]);
    c34[18] = fq_monomio(4, y2[7], y2[15], y1[19], y2[2]);
    c34[19] = fq_monomio(3, y4[2], y2[15], y1[21]);
    c34[20] = fq_monomio(4, y1[7], y2[15], y3[2], y1[21]);
    c34[21] = fq_monomio(4, y2[7], y2[15], y1[21], y2[2]);
    c34[22] = fq_monomio(4, y3[7], y2[15], y1[19], y1[2]);
    c34[23] = fq_monomio(4, y3[7], y2[15], y1[21], y1[2]);
    c34[24] = fq_monomio(3, y4[7], y2[11], y1[19]);
    c34[25] = fq_monomio(4, y4[7], y1[15], y1[19], y1[11]);
    c34[26] = fq_monomio(3, y4[7], y2[11], y1[21]);
    c34[27] = fq_monomio(4, y4[7], y1[15], y1[21], y1[11]);
    c34[28] = fq_monomio(3, y4[7], y2[15], y1[19]);
    c34[29] = fq_monomio(3, y4[7], y2[15], y1[21]);
    c34[30] = fq_monomio(2, y2[2], y1[11]);
    c34[31] = fq_monomio(3, y1[7], y1[11], y1[2]);
    c34[32] = fq_monomio(2, y2[2], y1[15]);
    c34[33] = fq_monomio(3, y1[7], y1[15], y1[2]);
    c34[34] = fq_monomio(2, y2[7], y1[11]);
    c34[35] = fq_monomio(2, y2[7], y1[15]);
    c56[0]  = fq_monomio(3, y4[3], y2[12], y1[23]);
    c56[1]  = fq_monomio(4, y1[8], y2[12], y3[3], y1[23]);
    c56[2]  = fq_monomio(4, y4[3], y1[16], y1[23], y1[12]);
    c56[3]  = fq_monomio(5, y1[8], y1[16], y3[3], y1[23], y1[12]);
    c56[4]  = fq_monomio(4, y2[8], y2[12], y1[23], y2[3]);
    c56[5]  = fq_monomio(5, y2[8], y1[16], y1[23], y2[3], y1[12]);
    c56[6]  = fq_monomio(3, y4[3], y2[12], y1[26]);
    c56[7]  = fq_monomio(4, y1[8], y2[12], y3[3], y1[26]);
    c56[8]  = fq_monomio(4, y4[3], y1[16], y1[26], y1[12]);
    c56[9]  = fq_monomio(5, y1[8], y1[16], y3[3], y1[26], y1[12]);
    c56[10] = fq_monomio(4, y2[8], y2[12], y1[26], y2[3]);
    c56[11] = fq_monomio(5, y2[8], y1[16], y1[26], y2[3], y1[12]);
    c56[12] = fq_monomio(4, y3[8], y2[12], y1[23], y1[3]);
    c56[13] = fq_monomio(5, y3[8], y1[16], y1[23], y1[12], y1[3]);
    c56[14] = fq_monomio(4, y3[8], y2[12], y1[26], y1[3]);
    c56[15] = fq_monomio(5, y3[8], y1[16], y1[26], y1[12], y1[3]);
    c56[16] = fq_monomio(3, y4[3], y2[16], y1[23]);
    c56[17] = fq_monomio(4, y1[8], y2[16], y3[3], y1[23]);
    c56[18] = fq_monomio(4, y2[8], y2[16], y1[23], y2[3]);
    c56[19] = fq_monomio(3, y4[3], y2[16], y1[26]);
    c56[20] = fq_monomio(4, y1[8], y2[16], y3[3], y1[26]);
    c56[21] = fq_monomio(4, y2[8], y2[16], y1[26], y2[3]);
    c56[22] = fq_monomio(4, y3[8], y2[16], y1[23], y1[3]);
    c56[23] = fq_monomio(4, y3[8], y2[16], y1[26], y1[3]);
    c56[24] = fq_monomio(3, y4[8], y2[12], y1[23]);
    c56[25] = fq_monomio(4, y4[8], y1[16], y1[23], y1[12]);
    c56[26] = fq_monomio(3, y4[8], y2[12], y1[26]);
    c56[27] = fq_monomio(4, y4[8], y1[16], y1[26], y1[12]);
    c56[28] = fq_monomio(3, y4[8], y2[16], y1[23]);
    c56[29] = fq_monomio(3, y4[8], y2[16], y1[26]);
    c78[0]  = fq_monomio(3, y4[4], y2[24], y1[13]);
    c78[1]  = fq_monomio(4, y1[9], y2[24], y3[4], y1[13]);
    c78[2]  = fq_monomio(4, y4[4], y1[27], y1[24], y1[13]);
    c78[3]  = fq_monomio(5, y1[9], y1[27], y3[4], y1[24], y1[13]);
    c78[4]  = fq_monomio(4, y2[9], y2[24], y2[4], y1[13]);
    c78[5]  = fq_monomio(5, y2[9], y1[27], y1[24], y2[4], y1[13]);
    c78[6]  = fq_monomio(3, y4[4], y2[27], y1[13]);
    c78[7]  = fq_monomio(4, y1[9], y2[27], y3[4], y1[13]);
    c78[8]  = fq_monomio(4, y2[9], y2[27], y2[4], y1[13]);
    c78[9]  = fq_monomio(4, y3[9], y2[24], y1[13], y1[4]);
    c78[10] = fq_monomio(5, y3[9], y1[27], y1[24], y1[13], y1[4]);
    c78[11] = fq_monomio(4, y3[9], y2[27], y1[13], y1[4]);
    c78[12] = fq_monomio(3, y4[4], y2[24], y1[17]);
    c78[13] = fq_monomio(4, y1[9], y2[24], y3[4], y1[17]);
    c78[14] = fq_monomio(4, y4[4], y1[27], y1[24], y1[17]);
    c78[15] = fq_monomio(5, y1[9], y1[27], y3[4], y1[24], y1[17]);
    c78[16] = fq_monomio(4, y2[9], y2[24], y2[4], y1[17]);
    c78[17] = fq_monomio(5, y2[9], y1[27], y1[24], y2[4], y1[17]);
    c78[18] = fq_monomio(3, y4[4], y2[27], y1[17]);
    c78[19] = fq_monomio(4, y1[9], y2[27], y3[4], y1[17]);
    c78[20] = fq_monomio(4, y2[9], y2[27], y2[4], y1[17]);
    c78[21] = fq_monomio(4, y3[9], y2[24], y1[17], y1[4]);
    c78[22] = fq_monomio(5, y3[9], y1[27], y1[24], y1[17], y1[4]);
    c78[23] = fq_monomio(4, y3[9], y2[27], y1[17], y1[4]);
    c78[24] = fq_monomio(3, y4[9], y2[24], y1[13]);
    c78[25] = fq_monomio(4, y4[9], y1[27], y1[24], y1[13]);
    c78[26] = fq_monomio(3, y4[9], y2[27], y1[13]);
    c78[27] = fq_monomio(3, y4[9], y2[24], y1[17]);
    c78[28] = fq_monomio(4, y4[9], y1[27], y1[24], y1[17]);
    c78[29] = fq_monomio(3, y4[9], y2[27], y1[17]);
    memset(ct, 0, 8*sizeof(fq_elem));
    for (i=0; i<78; i++)
    {
         ct[0] = fq_add(ct[0], fq_mul(c12[i], pk->p12[i][0]));
         ct[1] = fq_add(ct[1], fq_mul(c12[i], pk->p12[i][1]));
    }
    for (i=0; i<36; i++)
    {
         ct[2] = fq_add(ct[2], fq_mul(c34[i], pk->p34[i][0]));
         ct[3] = fq_add(ct[3], fq_mul(c34[i], pk->p34[i][1]));
    }
    for (i=0; i<30; i++)
    {
         ct[4] = fq_add(ct[4], fq_mul(c56[i], pk->p56[i][0]));
         ct[5] = fq_add(ct[5], fq_mul(c56[i], pk->p56[i][1]));
    }
    for (i=0; i<30; i++)
    {
         ct[6] = fq_add(ct[6], fq_mul(c78[i], pk->p78[i][0]));
         ct[7] = fq_add(ct[7], fq_mul(c78[i], pk->p78[i][1]));
    }
}

const unsigned char *pkey_parse(struct pkey_t *pk, const unsigned char *s)
{
    unsigned int i;
    memset(pk, 0, sizeof(struct pkey_t));
    for (i=0; i<78; i++)
        s = fq2_parse(pk->p12[i], s);
    for (i=0; i<36; i++)
        s = fq2_parse(pk->p34[i], s);
    for (i=0; i<30; i++)
        s = fq2_parse(pk->p56[i], s);
    for (i=0; i<30; i++)
        s = fq2_parse(pk->p78[i], s);
//    if ((pk->f[0] = *s++) >= 64) return NULL;
    if ((pk->f[1] = *s++) >= 64) return NULL;
    if ((pk->f[2] = *s++) >= 64) return NULL;
    if ((pk->f[3] = *s++) >= 64) return NULL;
    if ((pk->f[4] = *s++) >= 64) return NULL;
    if ((pk->f[5] = *s++) >= 64) return NULL;
//    if ((pk->f[6] = *s++) >= 64) return NULL;
    if ((pk->f[7] = *s++) >= 64) return NULL;
//    if ((pk->f[8] = *s++) >= 64) return NULL;
    if ((pk->f[9] = *s++) >= 64) return NULL;
//    if ((pk->f[10] = *s++) >= 64) return NULL;
    if ((pk->f[11] = *s++) >= 64) return NULL;
    if ((pk->f[12] = *s++) >= 64) return NULL;
//    if ((pk->f[13] = *s++) >= 64) return NULL;
    pk->f[0]  = (pk->f[1]  + delta) & 0x3f;
    pk->f[6]  = (pk->f[5]  + pk->f[2] - pk->f[1]) & 0x3f;
    pk->f[8]  = (pk->f[7]  + pk->f[4] - pk->f[3]) & 0x3f;
    pk->f[10] = (pk->f[9]  + pk->f[2] - pk->f[1]) & 0x3f;
    pk->f[13] = (pk->f[12] + pk->f[4] - pk->f[3]) & 0x3f;
    return s;
}

unsigned char *pkey_serialize(unsigned char *s, const struct pkey_t *pk)
{
    unsigned int i;
    for (i=0; i<78; i++)
        s = fq2_serialize(s, pk->p12[i]);
    for (i=0; i<36; i++)
        s = fq2_serialize(s, pk->p34[i]);
    for (i=0; i<30; i++)
        s = fq2_serialize(s, pk->p56[i]);
    for (i=0; i<30; i++)
        s = fq2_serialize(s, pk->p78[i]);
//    *s++ = pk->f[0];
    *s++ = pk->f[1];
    *s++ = pk->f[2];
    *s++ = pk->f[3];
    *s++ = pk->f[4];
    *s++ = pk->f[5];
//    *s++ = pk->f[6];
    *s++ = pk->f[7];
//    *s++ = pk->f[8];
    *s++ = pk->f[9];
//    *s++ = pk->f[10];
    *s++ = pk->f[11];
    *s++ = pk->f[12];
//    *s++ = pk->f[13];
    return s;    
}
