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

class ExpMean {
public:
  explicit ExpMean(double alpha,int numsamples)
  :alpha_(alpha),count_(numsamples),
   right_sum_(numsamples),inv_weight_(numsamples)
  {
    assert(numsamples>0 && alpha>0 && alpha<1.0);
    //calc inv_weight inplace
    //pass1: inv_weight holds right_weight[i]
    //right_weight[i]=alpha+alpha^2 + ... + alpha^(count_-1-i)
    inv_weight_[count_-1]=0.0;
    for (int i=count_-2;i>=0;--i) {
      inv_weight_[i]=alpha_*(1.0 + inv_weight_[i+1]);
    }
    //pass2: fold in left_weight[i]=alpha + ... + alpha^i
    double left_weight_=0;
    for (int i=0;i<count_;++i) {
      inv_weight_[i]=1.0/(left_weight_ + 1.0 + inv_weight_[i]);
      left_weight_ = alpha_ * (left_weight_ + 1.0);
    }
  }

  void Begin(const int32_t* buf,int bpn)
  {
    assert(bpn >= 0 && bpn < 32);
    buf_ = buf;
    pos_ = 0;
    left_sum_ = 0.0;

    const uint32_t plane = uint32_t{1} << bpn;
    known_mask_ = ~(plane - 1);
    high_mask_ = ~(plane | (plane - 1));

    right_sum_[count_-1]=0.0;
    for (int i = count_-2;i >= 0; --i)
      right_sum_[i] = alpha_*(High(i+1)+right_sum_[i+1]);
  }

  double Mean() const
  {
    assert(pos_>=0 && pos_<count_);
    return (left_sum_+High(pos_)+right_sum_[pos_])*inv_weight_[pos_];
  }

  void Advance()
  {
    assert(pos_<(count_-1));
    left_sum_ = alpha_ * (left_sum_ + Known(pos_));
    ++pos_;
  }

private:
  inline double Known(int i) const {
    return static_cast<uint32_t>(buf_[i]) & known_mask_;
  }

  inline double High(int i) const {
    return static_cast<uint32_t>(buf_[i]) & high_mask_;
  }

  double alpha_;
  int count_;
  std::vector<double> right_sum_, inv_weight_;
  const int32_t* buf_ = nullptr;
  int pos_ = 0;
  uint32_t known_mask_ = 0, high_mask_ = 0;
  double left_sum_ = 0;
};

class BitplaneCoder {
  const int cnt_upd_rate_p=150;
  const int cnt_upd_rate_sig=500;
  const int cnt_upd_rate_ref=200;
  const int cnt_upd_rate_sse=150;

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
    static int PredictLaplace(double mean,int bpn);
    int PredictRef();
    void UpdateRef(int bit);
    int PredictSig();
    void UpdateSig(int bit);
    int PredictSSE(int p1);
    void UpdateSSE(int bit);

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
