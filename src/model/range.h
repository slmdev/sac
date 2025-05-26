#ifndef RANGE_H
#define RANGE_H

#include "../common/bufio.h"
#include "model.h"
#include <functional>

class RangeCoderBase {
  public:
    explicit RangeCoderBase(BufIO &buf,int dec=0):buf(buf),decode(dec){};
    void SetDecode(){decode=1;};
    void SetEncode(){decode=0;};
  protected:
    BufIO &buf;
    int decode;
};

//#define SCALE_RANGE (((PSCALE-p1)*uint64_t(range)) >> PBITS) // 64 bit shift
#define SCALE_RANGE ((uint64_t(range)*((PSCALE-p1)<<(32-PBITS)))>>32)

#define RANGE_ENC_NORMALIZE  while ((low ^ (low+range))<TOP || (range<BOT && ((range= -(int)low & (BOT-1)),1))) buf.PutByte(low>>24),range<<=8,low<<=8;
#define RANGE_DEC_NORMALIZE  while ((low ^ (low+range))<TOP || (range<BOT && ((range= -(int)low & (BOT-1)),1))) (code<<=8)+=buf.GetByte(),range<<=8,low<<=8;

// Carryless RangeCoder
// derived from Dimitry Subbotin (public domain)
class RangeCoder : public RangeCoderBase
{
  enum : uint32_t {NUM=4,TOP=0x01000000U,BOT=0x00010000U};
  public:
    using RangeCoderBase::RangeCoderBase;
    void Init();
    void Stop();
    void EncodeSymbol(uint32_t low,uint32_t freq,uint32_t tot);
    void DecodeSymbol(uint32_t low,uint32_t freq);
    void EncodeBitOne(uint32_t p1,int bit);
    int  DecodeBitOne(uint32_t p1);
    uint32_t DecProb(uint32_t totfreq);
  protected:
    uint32_t low,range,code;
};


// Binary RangeCoder with Carry and 64-bit low
// derived from rc_v3 by Eugene Shelwien
class RangeCoderSH : public RangeCoderBase {
  enum : uint32_t { NUM=4,TOP=0x01000000U,Thres=0xFF000000U};
  public:
    using RangeCoderBase::RangeCoderBase;
    void Init();
    void Stop();
    void EncodeBitOne(uint32_t p1,int bit);
    int  DecodeBitOne(uint32_t p1);

    std::function<void(uint32_t,int)> encode_p1 = [this](uint32_t p1,int bit) {return EncodeBitOne(p1,bit);}; // stupid C++
    std::function<int(uint32_t)> decode_p1 = [this](uint32_t p1) {return DecodeBitOne(p1);};
  protected:
    void ShiftLow();
    uint32_t range,code,FFNum,Cache;
    uint64_t lowc;
};

#endif // RANGE_H
