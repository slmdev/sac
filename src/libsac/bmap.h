#ifndef BMAP_H
#define BMAP_H

#include "../global.h"
#include "../model/range.h"
#include "../model/counter.h"
#include "../model/mixer.h"
#include "../model/sse.h"

//encode binary maps from sparse-pcm model using context modeling

struct BIntMap {
  BIntMap() = default;
  BIntMap(span_u8 buf,int32_t minval,int32_t maxval)
  :buf(buf),minval(minval),maxval(maxval)
  {}
  inline int32_t idx2val(int32_t idx) const {return idx+minval;};
  inline int32_t val2idx(int32_t val) const {return val-minval;};
  span_u8 buf{};
  int32_t minval=0;
  int32_t maxval=0;
};

struct RunStats {
  int rb,rt,bhist;
  RunStats()
  :rb(0),rt(0),bhist(0){}
  void Update(int bit)
  {
    if (bit==rb) {if (rt<15) rt++;}
    else {rb=bit;rt=0;}

    bhist = ((bhist<<1)|bit)&31;
  }
};

class BMap {
  static constexpr int WCNT=750;
  static constexpr int WSSE=500;
  static constexpr int WMIX=750;
  public:
    BMap(RangeCoderSH &rc);
    void Encode(const BIntMap &bmap,const BIntMap &bmap_ref);
    void Decode(BIntMap &bmap,const BIntMap &bmap_ref);
  private:
    int32_t MirrorCtx(const BIntMap &bmap,int32_t idx) const;
    int32_t RefCtx(const BIntMap &bmap_ref,int32_t val) const;
    void Update(int bit);
    int32_t Predict();
    RangeCoderSH &rc;
    LinearCounter16 c0[32];
    LinearCounter16 c1[32];
    LinearCounter16 c2[32];
    LinearCounter16 c3[32];
    NMixLogistic pmix;
    SSENL<32> sse;
    int ctx0,ctx1,ctx2,ctx3;
    std::vector<int32_t>pv; //scratch-pad
};

#endif

