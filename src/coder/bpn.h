#ifndef BPN_H
#define BPN_H

#include <functional>
#include <bitset>
#include "pmodel.h"
#include "pstats.h"
#include "../model/counter.h"
#include "../model/sse.h"
#include "../model/mixer.h"
#include "../common/utils.h"
#include "apm2.h"

//#define h1y(v,k) (((v)>>k)^(v))
//#define h2y(v,k) (((v)*2654435761UL)>>(k))

template <int NBins>
class RelativeMagnitudeQ
{
  static_assert(NBins>=1);
  static constexpr int NB=NBins-1;
  public:
    RelativeMagnitudeQ(double growth)
    {
      assert(growth>1.0);
      double threshold=1.0/(growth*growth);
      for (double&t:thres) {
        t=threshold;
        threshold*=growth;
      }
    }
    int operator()(double r) const {
      assert(r>=0);
      int bin=0;
      while (bin<NB && r>=thres[bin])
        ++bin;
      return bin;
    }
  private:
    std::array<double,NB> thres;
};

class WeightedAgeQ
{
  public:
    WeightedAgeQ(int NBins,int num_ctx,double knee,double growth)
    :q_logbin(num_ctx)
    {
      for (int i=0;i<num_ctx;i++)
        q_logbin[i]=MathUtils::LogBin(i,knee,growth,NBins);
    }
    int operator()(int k) const {
      assert(k>=0 && k<(int)q_logbin.size());
      return q_logbin[k];
    }
  private:
    std::vector <int>q_logbin;
};

using EncodeP1 = std::function<void(uint32_t,int)>;
using DecodeP1 = std::function<int(uint32_t)>;

class BitplaneCoder {
  struct State {
    int pl_prob=0,logit0=0,logit0_q=0,logit1_q=0,logit2_q=0;
    int pl_scale=0;
    int sigst[17];
    int age_sum=0;
    int age_ctx=0;
    int rm_ctx=0;
    int ema_ctx1=0,ema_ctx2=0;
    int emb_ctx1=0,emb_ctx2=0;
    int nrun=0;
    uint32_t bhist=0;
  } s;

  static constexpr int limit_p=150;
  static constexpr int limit_sig=350;
  static constexpr int limit_ref=150;
  static constexpr int limit_sse_ref=500;
  static constexpr int limit_sse_sig=150;

  static constexpr int mix_upd_rate_ref =BM::MixerRate(0.0035);
  static constexpr int mix_upd_rate_sig =BM::MixerRate(0.0025);
  static constexpr int mix_upd_rate_msse=BM::MixerRate(0.0025);
  public:
    BitplaneCoder(int maxbpn,int numsamples);
    void Encode(EncodeP1 encode_p1,int32_t *buf);
    void Decode(DecodeP1 decode_p1,int32_t *buf);
  private:
    template <bool IsDecode,class Coder>
      void Run(Coder &&code,int32_t *buf);

    void GetState(const std::array<double,3> &means,uint32_t plane,int sample,State &s) const;

    inline int SigAge(int msb) const {
      return msb - (bpn+1);
    }

    int PredictRef();
    void UpdateRef(int bit);
    int PredictSig();
    void UpdateSig(int bit);
    int PredictSSERef(int p1);
    void UpdateSSERef(int bit);

    int PredictSSESig(int p1);
    void UpdateSSESig(int bit);

    APM2 <> apm_sig,apm_ref;

    std::vector<LinearCounterLimit> csig0,csig1,csig2,cref0,cref1,cref2;
    LinearCounterLimit pl;
    std::vector<LinearCounterLimit> chist;
    LinearCounterLimit *phist;
    std::vector <LogMixer>lmixref,lmixsig;

    std::vector<ScaleModel> sm;
    SigStats sigstats;
    RelativeMagnitudeQ<8> rm_q;
    WeightedAgeQ age_q;

    LinearCounterLimit *pc1,*pc2,*pc3,*pc4;
    LogMixer *plmix;
    int *pabuf,sample;
    //0=insignificant, else first significant plane+1
    std::vector <int>msb;
    uint32_t bmask[32];
    int maxbpn,bpn,numsamples;
    SigStats::St st;
};

#endif
