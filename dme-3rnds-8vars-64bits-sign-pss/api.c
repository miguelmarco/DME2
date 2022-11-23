#include <string.h>
#include "rng.h"
#include "dme.h"
#include "api.h"
#include "sha256.h"

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk)
{
    struct skey_t skey;
    struct pkey_t pkey;
    generate_skey(&skey);
    generate_pkey(&pkey, &skey);
    skey_serialize(sk, &skey);
    pkey_serialize(pk, &pkey);
    return 0;
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen,
    const unsigned char *m, unsigned long long mlen,
    const unsigned char *sk)
{
    unsigned int i;
    fq_elem mhash[8], signature[8];
    unsigned char r[16], w[32], g[32];
    struct skey_t skey;
    struct sha2_256_ctx ctx;
    if (!skey_parse(&skey, (unsigned char *) sk))
        return -1;
    randombytes(r, 16);
    sha2_256_init_ctx(&ctx);
    sha2_256_process_bytes(m, mlen, &ctx);
    sha2_256_process_bytes(r, 16,   &ctx);
    sha2_256_finish_ctx(&ctx, w);
//    if (*smlen < mlen + 64)
//        return -1;
    *smlen = mlen + 64;
    memcpy(sm, m, mlen);
    sha2_256_init_ctx(&ctx);
    sha2_256_process_bytes(w, 32, &ctx);
    sha2_256_finish_ctx(&ctx, g);
    for (i=0; i<16; i++)
        g[i] ^= r[i];
    memcpy(&mhash[0], w, 32);
    memcpy(&mhash[4], g, 32);
    decrypt_with_skey(signature, &skey, mhash);
    for (i=0; i<8; i++)
        fq_serialize(&sm[mlen+8*i], signature[i]);
    return 0;
}

int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
    const unsigned char *sm, unsigned long long smlen,
    const unsigned char *pk)
{
    unsigned int i;
    struct pkey_t pkey;
    struct sha2_256_ctx ctx;
    fq_elem signature[8], hash2[8];
    unsigned char g[32], w[32], r[16], g2[32], w2[32];
    if (!pkey_parse(&pkey, (unsigned char *) pk))
        return -1;
//    if (*mlen + 64 < smlen)
//        return -1;
    *mlen = smlen - 64;
    memcpy(m, sm, *mlen);
    for (i=0; i<8; i++)
        if (!fq_parse(&signature[i], &sm[*mlen + 8*i]))
            return -1;
    encrypt_with_pkey(hash2, &pkey, signature);
    memcpy(w, &hash2[0], 32);
    memcpy(g, &hash2[4], 32);
    sha2_256_init_ctx(&ctx);
    sha2_256_process_bytes(w, 32, &ctx);
    sha2_256_finish_ctx(&ctx, g2);
    for (i=0; i<16; i++)
        r[i] = g2[i] ^ g[i];
    if (memcmp(&g[16], &g2[16], 16))
        return -1;
    sha2_256_init_ctx(&ctx);
    sha2_256_process_bytes(m, *mlen, &ctx);
    sha2_256_process_bytes(r, 16, &ctx);
    sha2_256_finish_ctx(&ctx, w2);
    if (memcmp(w, w2, 32))
        return -1;
    return 0;
}
