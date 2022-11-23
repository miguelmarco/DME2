#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "api.h"
#include "rng.h"

#define TRIALS 1000
#define MAXBYTES 200

void hex(char *s, void *x, unsigned int n)
{
    unsigned int i;
    printf("%s", s);
    for (i=0; i<n; i++)
        printf("%.2x", ((unsigned char *)x)[i]);
    printf("\n");
}

int main(void)
{
    int i;
    clock_t t_gkey, t_sign, t_open, t1, t2;
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char m[MAXBYTES], m2[MAXBYTES];
    unsigned char sm[MAXBYTES + CRYPTO_BYTES];
    unsigned long long mlen, m2len, smlen;
    t_gkey = 0;
    t_sign = 0;
    t_open = 0;
    for (i=0; i<TRIALS; i++)
    {
        printf("--- TRIAL %u/%u ---\n", i+1, TRIALS);
        
        t1 = clock();
        if (crypto_sign_keypair(pk, sk))
        {
            printf("error: crypto_sign_keypair()\n");
            return -1;
        }
        t2 = clock();
        t_gkey += t2-t1;
        
        hex("sk = ", sk, CRYPTO_SECRETKEYBYTES);
        hex("pk = ", pk, CRYPTO_PUBLICKEYBYTES);

        randombytes(m, mlen = MAXBYTES);
        smlen = MAXBYTES + CRYPTO_BYTES;
        
        t1 = clock();
        if (crypto_sign(sm, &smlen, m, mlen, sk))
        {
            printf("error: crypto_sign()\n");
            return -1;
        }
        t2 = clock();
        t_sign += t2-t1;

        hex("m = ", m, mlen);
        hex("sm = ", sm, smlen);
        
        m2len = MAXBYTES;
        
        t1 = clock();
        if (crypto_sign_open(m2, &m2len, sm, smlen, pk))
        {
            printf("error: crypto_sign_open()\n");
            return -1;
        }
        t2 = clock();
        t_open += t2-t1;
        
        if (mlen != m2len || memcmp(m, m2, m2len))
        {
            printf("error: m != m2\n");
            hex("m = ", m2, m2len);
            return -1;
        }
        
        printf("\n");
    }
    
    printf("--- TIMINGS ---\n");
    printf("t_gkey: %.3f [usec]\n", t_gkey*1e6/CLOCKS_PER_SEC/TRIALS);
    printf("t_sign: %.3f [usec]\n", t_sign*1e6/CLOCKS_PER_SEC/TRIALS);
    printf("t_open: %.3f [usec]\n", t_open*1e6/CLOCKS_PER_SEC/TRIALS);
    
    return 0;
}
