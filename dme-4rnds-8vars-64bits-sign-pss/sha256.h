#ifndef SHA256_H
#define SHA256_H

#include <stdint.h>

uint32_t rrot32(uint32_t x, uint32_t c);
uint32_t get32be(const unsigned char *x);
void put32be(unsigned char *x, uint32_t r);
void put64be(unsigned char *x, uint64_t r);

struct sha2_256_ctx
{
  uint32_t h[8];
  unsigned char buf[64];
  size_t pending;
  uint64_t length;
};

void sha2_256_init_ctx(struct sha2_256_ctx *ctx);
void sha2_256_process_bytes(const unsigned char *buf, size_t len, struct sha2_256_ctx *ctx);
void sha2_256_finish_ctx(struct sha2_256_ctx *ctx, unsigned char *result);

#endif
