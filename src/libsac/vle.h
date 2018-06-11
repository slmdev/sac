#ifndef VLE_H
#define VLE_H

#include "../model/range.h"
#include "../model/counter.h"
#include "../model/sse.h"
#include "../model/mixer.h"
#include "../common/utils.h"

//#define h1y(v,k) (((v)>>k)^(v))
//#define h2y(v,k) (((v)*2654435761UL)>>(k))

class BitplaneCoder {
  const int cnt_upd_rate=350;
  const int mix_upd_rate=800;
  const int cntsse_upd_rate=1000;
  const int mixsse_upd_rate=300;
  public:
    BitplaneCoder(RangeCoderSH &rc,int maxbpn,int numsamples);
    void Encode(int32_t *abuf);
    void Decode(int32_t *buf);
  private:
    void GetSigState(int i); // get actual significance state
    int PredictRef();
    void UpdateRef(int bit);
    int PredictSig();
    void UpdateSig(int bit);
    int PredictSSE(int p1);
    void UpdateSSE(int bit);
    RangeCoderSH &rc;
    //int gc0[32],gc1[32];
    LinearCounter16 bpn1[1<<10],bpn2[1<<10];
    std::vector <NMixLogistic>lmixref,lmixsig;
    std::vector <NMixLogistic> finalmix;
    //SSENL<32> sse[2048];
    //SSENL<32> *psse;
    SSE<4> sse[2048];
    SSE<4> *psse;
    LinearCounter16 *pc1,*pc2,*pc3,*pc4,*pc5;
    NMixLogistic *plmix;
    int *pabuf,sample;
    std::vector <int>msb;
    int sigst[9];
    int maxbpn,bpn,numsamples;
    uint32_t bctx,state,sighist;
};

class Golomb {
  public:
    Golomb (RangeCoderSH &rc)
    :msum(0.8,1<<15),rc(rc)
    {
      lastl=0;
    }
    void Encode(int val)
    {
      if (val<0) val=2*(-val);
      else if (val>0) val=(2*val)-1;

      int m=std::max(static_cast<int>(msum.Get()),1);
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

      int m=std::max(static_cast<int>(msum.Get()),1);
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

      int m=std::max(static_cast<int>(msum.Get()),1);
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

      int m=std::max(static_cast<int>(msum.Get()),1);

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
