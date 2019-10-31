#ifndef MAP_H
#define MAP_H

#include "../model/range.h"
#include "../model/counter.h"
#include "../model/mixer.h"
#include "../model/sse.h"
#include <vector>

class MapEncoder {
  const int cnt_upd_rate=500;
  const int cntsse_upd_rate=300;
  const int mix_upd_rate=1000;
  const int mixsse_upd_rate=500;
  public:
    MapEncoder(RangeCoderSH &rc,std::vector <bool>&usedl,std::vector <bool>&usedh);
    void Encode();
    void Decode();
  private:
    int PredictLow(int i);
    int PredictHigh(int i);
    void Update(int bit);
    int PredictSSE(int p1,int ctx);
    void UpdateSSE(int bit,int ctx);
    RangeCoderSH &rc;
    LinearCounter16 cnt[24];
    LinearCounter16 cctx[256];
    LinearCounter16 *pc1,*pc2,*pc3,*pc4,*px;
    std::vector <NMixLogistic> mixl,mixh;
    NMixLogistic finalmix;
    NMixLogistic *mix;
    SSENL<32> sse[32];
    std::vector <bool>&ul,&uh;
    int lb;
};

class Remap {
  public:
    Remap();
    void Reset();
    double Compare(const Remap &cmap);
    void Analyse(int32_t *src,int numsamples);
    bool isUsed(int val);
    int32_t Map2(int32_t pred);
    int32_t Map(int32_t pred,int32_t err);
    int32_t Unmap(int32_t pred,int32_t merr);
    int scale,vmin,vmax;
    std::vector <bool>usedl,usedh;
    std::vector<int32_t> mapl,maph;
};


#endif // MAP_H
