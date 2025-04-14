#ifndef VLE_H
#define VLE_H

#include "../model/range.h"
#include "../model/counter.h"
#include "../model/sse.h"
#include "../model/mixer.h"
#include "../common/utils.h"
#include <functional>

//#define h1y(v,k) (((v)>>k)^(v))
//#define h2y(v,k) (((v)*2654435761UL)>>(k))

class StaticLaplaceModel {
  public:
    StaticLaplaceModel(int maxbpn)
    :pr((1<<maxbpn),std::vector<int>(32))
    {
      for (int sum=0;sum<(1<<maxbpn);sum++) {
        for (int bpn=0;bpn<32;bpn++) {
          double pd=0.;
          if (sum>0) {
            double theta=exp(-1.0/static_cast<double>(sum));
            pd=1.0-1.0/(1+pow(theta,1<<bpn));
          }
          int pi=std::clamp((int)round(pd*PSCALE),1,PSCALEm);
          pr[sum][bpn]=pi;
        }
      }
    }
    int Predict(int avg,int bpn)
    {
      return pr[avg][bpn];
    }
  private:
    std::vector<std::vector<int>> pr;
};

using EncodeP1 = std::function<void(uint32_t,int)>;
using DecodeP1 = std::function<int(uint32_t)>;

class BitplaneCoder {
  const int cnt_upd_rate_p=150;
  const int cnt_upd_rate_sig=300;
  const int cnt_upd_rate_ref=150;
  const int mix_upd_rate_ref=800;
  const int mix_upd_rate_sig=700;
  const int cntsse_upd_rate=250;
  const int mixsse_upd_rate=250;
  public:
    BitplaneCoder(int maxbpn,int numsamples);
    void Encode(EncodeP1 encode_p1,int32_t *abuf);
    void Decode(DecodeP1 decode_p1,int32_t *buf);
  private:
    void CountSig(int n,int &n1,int &n2);
    void GetSigState(int i); // get actual significance state
    int PredictLaplace(uint32_t avg_sum);
    int PredictRef();
    void UpdateRef(int bit);
    int PredictSig();
    void UpdateSig(int bit);
    int PredictSSE(int p1);
    void UpdateSSE(int bit);
    uint32_t GetAvgSum(int n);

    std::vector<LinearCounterLimit> csig0,csig1,csig2,csig3,cref0,cref1,cref2,cref3;
    std::vector<LinearCounterLimit>p_laplace;
    std::vector <NMixLogistic>lmixref,lmixsig;
    NMixLogistic ssemix;

    SSENL<15> sse[1<<12];
    SSENL<15> *psse1,*psse2;
    LinearCounterLimit *pc1,*pc2,*pc3,*pc4;
    LinearCounterLimit *pl;
    NMixLogistic *plmix;
    int *pabuf,sample;
    std::vector <int>msb;
    //int n_laplace;
    //std::vector <double>weights_laplace;
    int sigst[17];
    uint32_t bmask[32];
    int maxbpn,bpn,numsamples,nrun,pestimate;
    uint32_t state;
    StaticLaplaceModel lm;
};

class Golomb {
  public:
    Golomb (RangeCoderSH &rc)
    :msum(0.98,1<<15),rc(rc)
    {
      lastl=0;
    }
    void Encode(int val)
    {
      if (val<0) val=2*(-val);
      else if (val>0) val=(2*val)-1;

      int m=(std::max)(static_cast<int>(msum.sum),1);
      int q=val/m;
      int r=val-q*m;

      //for (int i=0;i<q;i++) rc.EncodeBitOne(PSCALEh,1); // encode exponent unary
      //rc.EncodeBitOne(PSCALEh,0);

      int ctx=1;
      for (int i=7;i>=0;i--) {
        int bit=(q>>i)&1;
        rc.EncodeBitOne(cnt[ctx].p1,bit);
        cnt[ctx].update(bit,250);

        ctx+=ctx+bit;
      }

      /*int ctx=0;
      for (int i=0;i<q;i++) {
        int pctx=lastl+(ctx<<1);
        rc.EncodeBitOne(cnt[pctx].p1,1);
        cnt[pctx].update(1,128);
        ctx++;
        if (ctx>1) ctx=1;
      }
      int pctx=lastl+(ctx<<1);
      rc.EncodeBitOne(cnt[pctx].p1,0);
      cnt[pctx].update(0,128);

      if (q>0) lastl=1;
      else lastl=0;*/

      if (m>1)
      {
        int b=ceil(log(m)/log(2));
        int t=(1<<b)-m;
        if (r < t) {
          for (int i=b-2;i>=0;i--) rc.EncodeBitOne(PSCALEh,((r>>i)&1));
        } else {
          for (int i=b-1;i>=0;i--) rc.EncodeBitOne(PSCALEh,(((r+t)>>i)&1));
        }
      }

      msum.Update(val);
    }
    int Decode() {
      int q=0;
      while (rc.DecodeBitOne(PSCALEh)!=0) q++;

      int m=(std::max)(static_cast<int>(msum.sum),1);
      int r=0;

      if (m>1)
      {
        int b=ceil(log(m)/log(2));
        int t=(1<<b)-m;
        for (int i=b-2;i>=0;i--) r=(r<<1)+rc.DecodeBitOne(PSCALEh);
        if (r>=t) r=((r<<1)+rc.DecodeBitOne(PSCALEh))-t;
      }

      int val=m*q+r;
      msum.Update(val);

      if (val) {
        if (val&1) val=((val+1)>>1);
        else val=-(val>>1);
      }
      return val;
    }
    RunExp msum;
  private:
    RangeCoderSH &rc;
    LinearCounter16 cnt[512];
    int lastl;
};

class GolombRC {
  public:
    GolombRC (RangeCoder &rc)
    :msum(0.8,1<<15),rc(rc)
    {
    }
    void Encode(int val)
    {
      if (val<0) val=2*(-val);
      else if (val>0) val=(2*val)-1;

      int m=(std::max)(static_cast<int>(msum.sum),1);
      int q=val/m;
      int r=val-q*m;

      for (int i=0;i<q;i++) rc.EncodeBitOne(PSCALEh,1); // encode exponent unary
      rc.EncodeBitOne(PSCALEh,0);

      rc.EncodeSymbol(r,1,m);

      msum.Update(val);
    }
    int Decode() {
      int q=0;
      while (rc.DecodeBitOne(PSCALEh)!=0) q++;

      int m=(std::max)(static_cast<int>(msum.sum),1);

      int r=rc.DecProb(m);
      rc.DecodeSymbol(r,1);

      int val=m*q+r;
      msum.Update(val);

      if (val) {
        if (val&1) val=((val+1)>>1);
        else val=-(val>>1);
      }
      return val;
    }
    RunExp msum;
  private:
    RangeCoder &rc;
};


#endif
