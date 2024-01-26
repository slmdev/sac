#ifndef MD5_H
#define MD5_H

#include <cstdint>

namespace MD5 {
  typedef struct{
    uint64_t size;        // Size of input in bytes
    uint32_t buffer[4];   // Current accumulation of hash
    uint8_t input[64];    // Input to be used in the next step
    uint8_t digest[16];   // Result of algorithm
  } MD5Context;

#define A 0x67452301
#define B 0xefcdab89
#define C 0x98badcfe
#define D 0x10325476


/*
 * Bit-manipulation functions defined by the MD5 algorithm
 */
#define F(X, Y, Z) ((X & Y) | (~X & Z))
#define G(X, Y, Z) ((X & Z) | (Y & ~Z))
#define H(X, Y, Z) (X ^ Y ^ Z)
#define I(X, Y, Z) (Y ^ (X | ~Z))

/*
 * Rotates a 32-bit word left by n bits
 */
  uint32_t rotateLeft(uint32_t x, uint32_t n);
  void Init(MD5Context *ctx);
  void Update(MD5Context *ctx, uint8_t *input, size_t input_len);
  void Finalize(MD5Context *ctx);
  void Step(uint32_t *buffer, uint32_t *input);
};

#endif
