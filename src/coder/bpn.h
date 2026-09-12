#ifndef BPN_H
#define BPN_H

#include <functional>
#include "../model/range.h"
#include "../model/counter.h"
#include "../model/sse.h"
#include "../model/mixer.h"
#include "../common/utils.h"

//#define h1y(v,k) (((v)>>k)^(v))
//#define h2y(v,k) (((v)*2654435761UL)>>(k))

using EncodeP1 = std::function<void(uint32_t,int)>;
using DecodeP1 = std::function<int(uint32_t)>;

struct AvgState {
  int max_radius;
  int scount = 0;
  uint64_t sum=0;
  uint32_t avg=0;
  explicit AvgState(int radius)
  :max_radius(radius) {};
};

class BitplaneCoder {
  const int cnt_upd_rate_p=150;
  const int cnt_upd_rate_sig=300;
  const int cnt_upd_rate_ref=150;
  const int cnt_upd_rate_sse=250;

  const int mix_upd_rate_ref=BM::MixerRate(0.005);
  const int mix_upd_rate_sig=BM::MixerRate(0.005);
  const int mix_upd_rate_sse=BM::MixerRate(0.0015);
  public:
    BitplaneCoder(int maxbpn,int numsamples);
    void Encode(EncodeP1 encode_p1,int32_t *abuf);
    void Decode(DecodeP1 decode_p1,int32_t *buf);
  private:
    void CountSig(int n,int &n1,int &n2);
    void GetSigState(int i); // get actual significance state
    int PredictLaplace(const AvgState &ctx);
    int PredictRef();
    void UpdateRef(int bit);
    int PredictSig();
    void UpdateSig(int bit);
    int PredictSSE(int p1);
    void UpdateSSE(int bit);
    void UpdateAvgState(AvgState &state);

    std::vector<LinearCounterLimit> csig0,csig1,csig2,csig3,cref0,cref1,cref2,cref3;
    std::vector<LinearCounterLimit>p_laplace;
    std::vector <LogMixer>lmixref,lmixsig;
    LogMixer ssemix;

    SSENL<15> sse[1<<12];
    SSENL<15> *psse1,*psse2;
    LinearCounterLimit *pc1,*pc2,*pc3,*pc4;
    LinearCounterLimit *pl;
    LogMixer *plmix;
    int *pabuf,sample;
    std::vector <int>msb;
    int sigst[17];
    uint32_t bmask[32];
    int maxbpn,bpn,numsamples,nrun,pestimate;
    uint32_t state;
};

#endif
