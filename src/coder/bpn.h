#ifndef BPN_H
#define BPN_H

#include <functional>
#include "pmodel.h"
#include "pstats.h"
#include "../model/counter.h"
#include "../model/sse.h"
#include "../model/mixer.h"

//#define h1y(v,k) (((v)>>k)^(v))
//#define h2y(v,k) (((v)*2654435761UL)>>(k))

using EncodeP1 = std::function<void(uint32_t,int)>;
using DecodeP1 = std::function<int(uint32_t)>;

class BitplaneCoder {
  static constexpr int cnt_upd_rate_p=125;
  static constexpr int cnt_upd_rate_sig=500;
  static constexpr int cnt_upd_rate_ref=200;
  static constexpr int cnt_upd_rate_sse=150;

  static constexpr int mix_upd_rate_ref=BM::MixerRate(0.005);
  static constexpr int mix_upd_rate_sig=BM::MixerRate(0.0045);
  static constexpr int mix_upd_rate_sse=BM::MixerRate(0.0015);

  static constexpr int SIG_RADIUS=32;
  static constexpr int SIG_AGE=2;
  public:
    BitplaneCoder(int maxbpn,int numsamples);
    void Encode(EncodeP1 encode_p1,int32_t *buf);
    void Decode(DecodeP1 decode_p1,int32_t *buf);
  private:
    template <bool IsDecode,class Coder>
      void Run(Coder &&code,int32_t *buf);

    static inline int LaplaceLogitBin(int x, int nbits);
    static inline int LogBin(unsigned x,int nbins);
    static inline int LogBin(unsigned x,double knee,double growth,int nbins);
    std::tuple<int,int,int> CountSig(int n);
    void GetSigState(int i); // get actual significance state
    int PredictRef();
    void UpdateRef(int bit);
    int PredictSig();
    void UpdateSig(int bit);
    int PredictSSE(int p1);
    void UpdateSSE(int bit);

    std::vector<LinearCounterLimit> csig0,csig1,csig2,csig3,cref0,cref1,cref2,cref3;
    LinearCounterLimit pl;
    std::vector <LogMixer>lmixref,lmixsig;
    LogMixer ssemix;
    std::vector<ScaleModel> sm;
    SigStats sigstats;

    SSENL<16> sse1[1<<8];
    SSENL<16> sse2[1<<8];
    SSENL<16> *psse1,*psse2;
    LinearCounterLimit *pc1,*pc2,*pc3,*pc4;
    LogMixer *plmix;
    int *pabuf,sample;
    std::vector <int>msb;
    int sigst[17];
    uint32_t bmask[32];
    int maxbpn,bpn,numsamples,nrun;
    int pl_prob,pl_logit,pl_scale;
    uint32_t state;
};

#endif
