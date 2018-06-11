#include "vle.h"

BitplaneCoder::BitplaneCoder(RangeCoderSH &rc,int maxbpn,int numsamples)
:rc(rc),lmixref(32,NMixLogistic(3)),lmixsig(32,NMixLogistic(2)),
finalmix(32,NMixLogistic(2)),
msb(numsamples),maxbpn(maxbpn),numsamples(numsamples)
{
  bctx=state=0;
  bpn=0;
}
void BitplaneCoder::GetSigState(int i)
{
  sigst[0]=msb[i];
  sigst[1]=i>0?msb[i-1]:0;
  sigst[2]=i<numsamples-1?msb[i+1]:0;
  sigst[3]=i>1?msb[i-2]:0;
  sigst[4]=i<numsamples-2?msb[i+2]:0;
  sigst[5]=i>2?msb[i-3]:0;
  sigst[6]=i<numsamples-3?msb[i+3]:0;
  sigst[7]=i>3?msb[i-4]:0;
  sigst[8]=i<numsamples-4?msb[i+4]:0;
}

int BitplaneCoder::PredictRef()
{
  int mixctx=(state&7);

  int val=pabuf[sample];

  int lval=sample>0?pabuf[sample-1]:0;
  int lval2=sample>1?pabuf[sample-2]:0;
  int nval=sample<(numsamples-1)?pabuf[sample+1]:0;
  int nval2=sample<(numsamples-2)?pabuf[sample+2]:0;

  int lctx=0;
  int tval=val>>(bpn+1);
  int tval1=tval<<1;
  if ((lval>>bpn)>tval1) lctx+=(1<<0);
  if (nval>>(bpn+1)>tval) lctx+=(1<<1);
  if ((lval2>>bpn)>tval1) lctx+=(1<<2);
  if (nval2>>(bpn+1)>tval) lctx+=(1<<3);
      //if (lval>>bpn>lval2>>bpn) lctx+=(1<<4);

  pc1=&bpn2[0];
  pc2=&bpn1[1+msb[sample]];
  pc3=&bpn1[1+32+lctx];
  plmix=&lmixref[mixctx];

  int px=plmix->Predict({pc1->p1,pc2->p1,pc3->p1});
  return px;
}

void BitplaneCoder::UpdateRef(int bit)
{
  pc1->update(bit,cnt_upd_rate);
  pc2->update(bit,cnt_upd_rate);
  pc3->update(bit,cnt_upd_rate);
  plmix->Update(bit,mix_upd_rate);
  state=(state<<1)+0;
}

int BitplaneCoder::PredictSig()
{
  int mixctx=(state&7);
  int ctx1=0;
      //int ctx2=0;

  if (sigst[1]) ctx1+=(1<<0);
  if (sigst[2]) ctx1+=(1<<1);
  if (sigst[3]) ctx1+=(1<<2);
  if (sigst[4]) ctx1+=(1<<3);
  if (sigst[5]) ctx1+=(1<<4);
  if (sigst[6]) ctx1+=(1<<5);
  if (sigst[7]) ctx1+=(1<<6);
  if (sigst[8]) ctx1+=(1<<7);

  pc1=&bpn2[0];
  pc2=&bpn2[1+ctx1];

  plmix=&lmixsig[mixctx];
  return plmix->Predict({pc1->p1,pc2->p1});
}

void BitplaneCoder::UpdateSig(int bit)
{
  pc1->update(bit,cnt_upd_rate);
  pc2->update(bit,cnt_upd_rate);
  plmix->Update(bit,mix_upd_rate);
  state=(state<<1)+1;
  sighist=(sighist<<1)+bit;
}

int BitplaneCoder::PredictSSE(int p1)
{
  int ssectx=(bpn<<5)+msb[sample];
  psse=&sse[ssectx];
  return finalmix[0].Predict({psse->Predict(p1),p1});
}

void BitplaneCoder::UpdateSSE(int bit)
{
  psse->Update(bit,cntsse_upd_rate);
  finalmix[0].Update(bit,mixsse_upd_rate);
}

void BitplaneCoder::Encode(int32_t *abuf)
{
  pabuf=abuf;
  for (bpn=maxbpn;bpn>=0;bpn--)  {
    bctx=1;
    state=0;
    sighist=1;
    for (sample=0;sample<numsamples;sample++) {
      GetSigState(sample);
      int bit=(pabuf[sample]>>bpn)&1;
      if (sigst[0]) { // coef is significant, refine
        //rc.EncodeBitOne(PredictRef(),bit);
        rc.EncodeBitOne(PredictSSE(PredictRef()),bit);
        UpdateRef(bit);
        UpdateSSE(bit);
      } else { // coef is insignificant
        //rc.EncodeBitOne(PredictSig(),bit);
        rc.EncodeBitOne(PredictSSE(PredictSig()),bit);
        UpdateSig(bit);
        UpdateSSE(bit);
        if (bit) msb[sample]=bpn;
      }
      bctx=((bctx<<1)+bit)&15;
    }
  }
}

void BitplaneCoder::Decode(int32_t *buf)
{
  int bit;
  pabuf=buf;
  for (int i=0;i<numsamples;i++) buf[i]=0;
  for (bpn=maxbpn;bpn>=0;bpn--)  {
    bctx=1;
    state=0;
    for (sample=0;sample<numsamples;sample++) {
      GetSigState(sample);
      if (sigst[0]) { // coef is significant, refine
        bit=rc.DecodeBitOne(PredictSSE(PredictRef()));
        UpdateRef(bit);
        UpdateSSE(bit);
        if (bit) buf[sample]+=(1<<bpn);
       } else { // coef is insignificant
         bit=rc.DecodeBitOne(PredictSSE(PredictSig()));
         UpdateSig(bit);
         UpdateSSE(bit);
         if (bit) {
           buf[sample]+=(1<<bpn);
           msb[sample]=bpn;
          }
        }
    }
  }
  for (int i=0;i<numsamples;i++) buf[i]=MathUtils::U2S(buf[i]);
}
