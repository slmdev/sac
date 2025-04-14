#pragma once

#include <cstdint>
#include <stddef.h>

/*
 * Bit-manipulation functions defined by the MD5 algorithm
 */
inline uint32_t F(uint32_t &X, uint32_t &Y, uint32_t &Z) {return ((X & Y) | (~X & Z));}
inline uint32_t G(uint32_t &X, uint32_t &Y, uint32_t &Z) {return ((X & Z) | (Y & ~Z));}
inline uint32_t H(uint32_t &X, uint32_t &Y, uint32_t &Z) {return (X ^ Y ^ Z);}
inline uint32_t I(uint32_t &X, uint32_t &Y, uint32_t &Z) {return (Y ^ (X | ~Z));}

namespace MD5 {
  struct MD5Context{
    uint64_t size;        // Size of input in bytes
    uint32_t buffer[4];   // Current accumulation of hash
    uint8_t input[64];    // Input to be used in the next step
    uint8_t digest[16];   // Result of algorithm
  };
  
  const uint32_t MD5_A = 0x67452301;
  const uint32_t MD5_B = 0xefcdab89;
  const uint32_t MD5_C = 0x98badcfe;
  const uint32_t MD5_D = 0x10325476;
  
  /*
   * Rotates a 32-bit word left by n bits
   */
  uint32_t rotateLeft(uint32_t x, uint32_t n);
  void Init(MD5Context *ctx);
  void Update(MD5Context *ctx, uint8_t *input, size_t input_len);
  void Finalize(MD5Context *ctx);
  void Step(uint32_t *buffer, uint32_t *input);
};
