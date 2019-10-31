#include "range.h"

#define DO(n) for (uint32_t _=0;_<n;_++)

void RangeCoder::Init()
{
  low     = code  = 0;
  range   = uint32_t(-1);
  if (decode==1) DO(NUM) (code <<=8) += buf.GetByte();
}

void RangeCoder::Stop()
{
  if (decode==0) DO(NUM) buf.PutByte(low>>24),low<<=8;
}

void RangeCoder::EncodeSymbol(uint32_t cumfreq,uint32_t freq,uint32_t totfreq)
{
  low   += cumfreq * (range /= totfreq);
  range *= freq;
  RANGE_ENC_NORMALIZE
}

void RangeCoder::DecodeSymbol(uint32_t cumfreq,uint32_t freq)
{
  low   += cumfreq*range;
  range *= freq;
  RANGE_DEC_NORMALIZE
}

uint32_t RangeCoder::DecProb(uint32_t totfreq)
{
  uint32_t tmp=(code-low) / (range /= totfreq);
  return tmp;
}

void RangeCoder::EncodeBitOne(uint32_t p1,const int bit)
{
  const uint32_t rnew=SCALE_RANGE;
  bit ? low += rnew,range-=rnew:range=rnew;
  RANGE_ENC_NORMALIZE
}

int RangeCoder::DecodeBitOne(uint32_t p1)
{
  const uint32_t rnew=SCALE_RANGE;
  int bit=(code-low>=rnew);
  bit ? low += rnew,range-=rnew:range=rnew;
  RANGE_DEC_NORMALIZE
  return bit;
}

// binary rangecoder
void RangeCoderSH::Init()
{
  range = 0xFFFFFFFF;
  lowc = FFNum = Cache = code = 0;
  if(decode==1) DO(NUM+1) (code <<=8) += buf.GetByte();
}

void RangeCoderSH::Stop()
{
  if (decode==0) DO(NUM+1) ShiftLow();
}

void RangeCoderSH::EncodeBitOne(uint32_t p1,int bit)
{
  const uint32_t rnew=SCALE_RANGE;
  bit ? range-=rnew, lowc+=rnew : range=rnew;
  while(range<TOP) range<<=8,ShiftLow();
}

int RangeCoderSH::DecodeBitOne(uint32_t p1)
{
  const uint32_t rnew=SCALE_RANGE;
  int bit = (code>=rnew);
  bit ? range-=rnew, code-=rnew : range=rnew;
  while(range<TOP) range<<=8,(code<<=8)+=buf.GetByte();
  return bit;
}

void RangeCoderSH::ShiftLow()
{
  uint32_t Carry = uint32_t(lowc>>32), low = uint32_t(lowc);
  if( low<Thres || Carry )
  {
     buf.PutByte(Cache+Carry);
     for (;FFNum != 0;FFNum--) buf.PutByte(Carry-1);
     Cache = low>>24;
   } else FFNum++;
  lowc = (low<<8);
}
