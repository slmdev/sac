#include "utils.h"


namespace BitUtils
{
  uint32_t get32HL(const uint8_t *buf)
  {
    return((uint32_t)buf[3] + ((uint32_t)buf[2] << 8) +((uint32_t)buf[1] << 16) + ((uint32_t)buf[0] << 24));
  }
  uint16_t get16LH(const uint8_t *buf)
  {
    return((uint16_t)buf[0] + ((uint16_t)buf[1] << 8));
  }
  uint32_t get32LH(const uint8_t *buf)
  {
    return((uint32_t)buf[0] + ((uint32_t)buf[1] << 8) +((uint32_t)buf[2] << 16) + ((uint32_t)buf[3] << 24));
  }
  void put16LH(uint8_t *buf,uint16_t val)
  {
    buf[0] = val & 0xff;
    buf[1] = (val>>8) & 0xff;
  }
  void put32LH(uint8_t *buf,uint32_t val)
  {
    buf[0] = val & 0xff;
    buf[1] = (val>>8) & 0xff;
    buf[2] = (val>>16) & 0xff;
    buf[3] = (val>>24) & 0xff;
  }
  std::string U322Str(uint32_t val)
  {
    std::string s;
    for (int i=0;i<4;i++) {s+=(char)(val & 0xff);val>>=8;};
    return s;
  }
}
