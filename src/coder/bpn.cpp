#include "bpn.h"

BitplaneCoder::BitplaneCoder(int maxbpn,int numsamples)
:csig0(1<<20),csig1(1<<20),csig2(1<<20),csig3(1<<20),
cref0(1<<20),cref1(1<<20),cref2(1<<20),cref3(1<<20),
lmixref(256,LogMixer(5)),lmixsig(256,LogMixer(3)),
ssemix(2),
msb(numsamples),
maxbpn(maxbpn),numsamples(numsamples)
{
  state=0;
  bpn=0;
  nrun=0;
  pestimate=0;
  for (int i=0;i<32;i++)
    bmask[i]=~((uint32_t{1}<<i)-1);
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
  sigst[9]=i>4?msb[i-5]:0;
  sigst[10]=i<numsamples-5?msb[i+5]:0;
  sigst[11]=i>5?msb[i-6]:0;
  sigst[12]=i<numsamples-6?msb[i+6]:0;
  sigst[13]=i>6?msb[i-7]:0;
  sigst[14]=i<numsamples-7?msb[i+7]:0;
  sigst[15]=i>7?msb[i-8]:0;
  sigst[16]=i<numsamples-8?msb[i+8]:0;
}


int BitplaneCoder::PredictLaplace(double mean,int bpn)
{
  if (mean<0.1) return 1;
  const double plane=static_cast<double>(uint32_t{1}<<bpn);
  double t=std::exp(-plane/mean);
  double ps=(t/(1.0+t))*PSCALE;
  return std::clamp((int)std::round(ps),1,PSCALEm);
}

int BitplaneCoder::PredictRef()
{
  int val=pabuf[sample];

  int lval=sample>0?pabuf[sample-1]:0;
  int lval2=sample>1?pabuf[sample-2]:0;
  int nval=sample<(numsamples-1)?pabuf[sample+1]:0;
  int nval2=sample<(numsamples-2)?pabuf[sample+2]:0;

  int b0=(val>>(bpn+1));
  int b1=(lval>>(bpn));
  int b2=(nval>>(bpn+1));
  int b3=(lval2>>(bpn));
  int b4=(nval2>>(bpn+1));


  int c0=(b0<<1)<b1?1:0;
  int c1=(b0)<b2?1:0;
  int c2=(b0<<1)<b3?1:0;
  int c3=(b0)<(b4)?1:0;


  int x0=(val>>(bpn+1))<<1;
  int x1=(lval>>bpn);
  int x2=(nval>>(bpn+1))<<1;
  int x3=(lval2>>(bpn));
  int x4=(nval2>>(bpn+1))<<1;
  int xm=(x0+x1+x2+x3+x4)/5;

  int d0=x0>xm;
  //int d1=x1>xm;
  //int d2=x2>xm;

  int ctx1=(b0&31)+((b1&7)<<5)+((b3&1)<<8); //9-bits
  int ctx2=(c0+(c1<<1)+(c2<<2)+(c3<<3))+(d0<<4); //5-bits

  pc1=&cref0[msb[sample]-bpn];
  pc2=&cref1[ctx1];
  pc3=&cref2[ctx2];

  int ref_mixctx=((msb[sample])<<1)+d0;
  plmix=&lmixref[ref_mixctx];

  int px=plmix->Predict({pestimate,pl.p1,pc1->p1,pc2->p1,pc3->p1});

  return px;
}

void BitplaneCoder::UpdateRef(int bit)
{
  pl.update(bit,cnt_upd_rate_p);
  pc1->update(bit,cnt_upd_rate_ref);
  pc2->update(bit,cnt_upd_rate_ref);
  pc3->update(bit,cnt_upd_rate_ref);
  plmix->Update(bit,mix_upd_rate_ref);
  state=(state<<1)+0;
}

// count number of significant samples in neighborhood
std::tuple<int,int> BitplaneCoder::CountSig(int n)
{
  int n1=0,n2=0;

  auto update_scores=[&](int x)
  {
    if (x) {
      const int d=std::min(x-bpn,2);
      n1+=1;
      n2+=(d!=0);
    }
  };

  for (int i=1;i<=n;i++) {
    if (sample-i>=0)
      update_scores(msb[sample-i]);

    if (sample+i<numsamples)
      update_scores(msb[sample+i]);
  }
  return {n1,n2};
}

int BitplaneCoder::PredictSig()
{
  int ctx1=0;
  for (int i=0;i<16;i++)
    if (sigst[i+1]) ctx1+=1<<i;

  const auto [n1,n2] = CountSig(32);
  int ctx2=n2;

  pc1=&csig0[ctx1];
  pc2=&csig1[ctx2];

  //8-bit
  int sig_mixctx=(nrun<<4)+((n1>=3?3:n1)<<2)+((n2>0)<<1)+(n1-n2>1);
  plmix=&lmixsig[sig_mixctx];
  int p_mix=plmix->Predict({pl.p1,pc1->p1,pc2->p1});
  return p_mix;
}

void BitplaneCoder::UpdateSig(int bit)
{
  pl.update(bit,cnt_upd_rate_p);
  pc1->update(bit,cnt_upd_rate_sig);
  pc2->update(bit,cnt_upd_rate_sig);
  plmix->Update(bit,mix_upd_rate_sig);
  state=(state<<1)+1;
}

int BitplaneCoder::PredictSSE(int p1)
{
  int dist=sigst[0]?std::min(msb[sample]-bpn,7):0;
  int ctx1=((pestimate>>(PBITS-4))<<3)+(dist);
  int ctx2=(sigst[0]?1:0)+((sigst[1]?1:0)<<1)+((sigst[2]?1:0)<<2)+((sigst[3]?1:0)<<3)+((sigst[4]?1:0)<<4)+((sigst[5]?1:0)<<5)+((sigst[6]?1:0)<<6);
  psse1=&sse1[ctx1];
  psse2=&sse2[ctx2];
  int pr1=psse1->Predict(p1);
  int pr2=psse2->Predict(pr1);
  return ssemix.Predict({(pr1+pr2+1)>>1,p1});
}

void BitplaneCoder::UpdateSSE(int bit)
{
  psse1->Update(bit,cnt_upd_rate_sse);
  psse2->Update(bit,cnt_upd_rate_sse);
  ssemix.Update(bit,mix_upd_rate_sse);
}

void BitplaneCoder::Encode(EncodeP1 encode_p1,int32_t *abuf)
{
  pabuf=abuf;
  ExpMean em(0.95,numsamples);
  for (bpn=maxbpn;bpn>=0;bpn--)  {
    em.Begin(pabuf,bpn);
    state=0;
    nrun=0;
    pl.p1 = PredictLaplace(em.Mean(),bpn);
    for (sample=0;sample<numsamples;sample++) {
      pestimate=PredictLaplace(em.Mean(),bpn);
      GetSigState(sample);
      int bit=(pabuf[sample]>>bpn)&1;
      int p=0;
      if (sigst[0]) { // coef is significant, refine
        nrun=0;
        p=PredictSSE(PredictRef());
        encode_p1(p,bit);
        UpdateRef(bit);
        UpdateSSE(bit);
      } else { // coef is insignificant
        if (nrun<15) nrun++;
        p=PredictSSE(PredictSig());
        encode_p1(p,bit);
        UpdateSig(bit);
        UpdateSSE(bit);
        if (bit) msb[sample]=bpn;
      }
      if (sample<numsamples-1) em.Advance();
    }
  }
}

void BitplaneCoder::Decode(DecodeP1 decode_p1,int32_t *buf)
{
  int bit;
  pabuf=buf;
  ExpMean em(0.95,numsamples);
  for (int i=0;i<numsamples;i++) buf[i]=0;
  for (bpn=maxbpn;bpn>=0;bpn--)  {
    em.Begin(pabuf,bpn);
    state=0;
    nrun=0;
    pl.p1 = PredictLaplace(em.Mean(),bpn);
    for (sample=0;sample<numsamples;sample++) {
      pestimate=PredictLaplace(em.Mean(),bpn);
      GetSigState(sample);
      if (sigst[0]) { // coef is significant, refine
        nrun=0;
        bit=decode_p1(PredictSSE(PredictRef()));
        UpdateRef(bit);
        UpdateSSE(bit);
        if (bit) buf[sample]+=(1<<bpn);
      } else { // coef is insignificant
         if (nrun<15) nrun++;
         bit=decode_p1(PredictSSE(PredictSig()));
         UpdateSig(bit);
         UpdateSSE(bit);
         if (bit) {
           buf[sample]+=(1<<bpn);
           msb[sample]=bpn;
          }
      }
      if (sample<numsamples-1) em.Advance();
    }
  }
  for (int i=0;i<numsamples;i++) buf[i]=MathUtils::U2S(buf[i]);
}

