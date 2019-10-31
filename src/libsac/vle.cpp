#include "vle.h"

BitplaneCoder::BitplaneCoder(RangeCoderSH &rc,int maxbpn,int numsamples)
:rc(rc),
csig0(1<<20),csig1(1<<20),csig2(1<<20),csig3(1<<20),
cref0(1<<20),cref1(1<<20),cref2(1<<20),cref3(1<<20),
p_laplace(32),
lmixref(256,NMixLogistic(6)),lmixsig(256,NMixLogistic(3)),
ssemix(2),
msb(numsamples),maxbpn(maxbpn),numsamples(numsamples)
{
  state=0;
  bpn=0;
  nrun=0;
  double theta=0.99;
  for (int i=0;i<32;i++) {
    int p=std::min(std::max((int)round((1-1.0/(1+pow(theta,1<<i)))*PSCALE),1),PSCALEm);
    //std::cout << p << ' ';
    p_laplace[i].p1=p;
  }
  pestimate=0;
  for (int i=0;i<32;i++) {
    bmask[i]=~((1<<i)-1);
  }
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

static inline uint32_t ilog2(const uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}

int BitplaneCoder::PredictLaplace()
{
  int n=32;
  int nsum=0;
  int nidx=0;
  for (int k=sample-n;k<=sample+n;k++) {
    if (k>=0 && k<numsamples) {
      int val=pabuf[k];
      val&=k<sample?bmask[bpn]:bmask[bpn+1];
      nsum+=val;
      nidx++;
    }
  }
  double zm=nsum>0?double(nidx)/nsum:0.;
  double theta=exp(-zm);
  double p_l=1.0-1.0/(1+pow(theta,1<<bpn));
  return std::min(std::max((int)round(p_l*PSCALE),1),PSCALEm);
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
  int d1=x1>xm;
  //int d2=x2>xm;

  int ctx1=(b0&15)+((b1&15)<<4)+((b2&15)<<8);
  int ctx2=(c0+(c1<<1)+(c2<<2)+(c3<<3))+(d0<<4)+(d1<<5);
  int ctx3=(sigst[1]+sigst[2]+sigst[3]+sigst[4]+sigst[5]+sigst[6]+sigst[7]+sigst[8]);

  pl=&p_laplace[bpn];
  pc1=&cref0[msb[sample]];
  pc2=&cref1[ctx1&255];
  pc3=&cref2[ctx2];
  pc4=&cref3[ctx3];

  //int mixctx=((state&15)<<1)+d0;

  int pctx=((((pestimate>>12)<<1)+d0)<<1)+(b0&1);

  plmix=&lmixref[pctx];
  int px=plmix->Predict({pestimate,pl->p1,pc1->p1,pc2->p1,pc3->p1,pc4->p1});
  return px;
}

void BitplaneCoder::UpdateRef(int bit)
{
  pl->update(bit,cnt_upd_rate_p);
  pc1->update(bit,cnt_upd_rate_ref);
  pc2->update(bit,cnt_upd_rate_ref);
  pc3->update(bit,cnt_upd_rate_ref);
  pc4->update(bit,cnt_upd_rate_ref);
  plmix->Update(bit,mix_upd_rate_ref);
  state=(state<<1)+0;
}

// count number of significant samples in neighborhood
void BitplaneCoder::CountSig(int n,int &n1,int &n2)
{
  n1=n2=0;
  for (int i=1;i<=n;i++) {
    if (sample-i>=0) {
       if (msb[sample-i]) n1+=1;
       if (msb[sample-i]>bpn) n2+=1;
    }
    if (sample+i<numsamples-1) {
       if (msb[sample+i]) n1+=1;
       if (msb[sample+i]>bpn) n2+=1;
    }
  }
}

int BitplaneCoder::PredictSig()
{
  int ctx1=0;
  for (int i=0;i<16;i++)
    if (sigst[i+1]) ctx1+=1<<i;

  int n1,n2;
  CountSig(32,n1,n2);
  int ctx2=n2;

  pl=&p_laplace[bpn];
  pc1=&csig0[ctx1];
  pc2=&csig1[ctx2];

  int mixctx=((state&15)<<3)+((n1>=3?3:n1)<<1)+(n2>0?1:0);
  plmix=&lmixsig[mixctx];
  int p_mix=plmix->Predict({pl->p1,pc1->p1,pc2->p1});
  return p_mix;
}

void BitplaneCoder::UpdateSig(int bit)
{
  pl->update(bit,cnt_upd_rate_p);
  pc1->update(bit,cnt_upd_rate_sig);
  pc2->update(bit,cnt_upd_rate_sig);
  plmix->Update(bit,mix_upd_rate_sig);
  state=(state<<1)+1;
}

int BitplaneCoder::PredictSSE(int p1)
{
  //((pestimate>>11)<<1)+
  psse1=&sse[((pestimate>>11)<<1)+(sigst[0]?1:0)];
  psse2=&sse[32+(sigst[0]?1:0)+((sigst[1]?1:0)<<1)+((sigst[2]?1:0)<<2)+((sigst[3]?1:0)<<3)+((sigst[4]?1:0)<<4)+((sigst[5]?1:0)<<5)+((sigst[6]?1:0)<<6)];
  int pr1=psse1->Predict(p1);
  int pr2=psse2->Predict(pr1);

  return ssemix.Predict({(pr1+pr2+1)>>1,p1});
}

void BitplaneCoder::UpdateSSE(int bit)
{
  psse1->Update(bit,cntsse_upd_rate,true);
  psse2->Update(bit,cntsse_upd_rate,true);
  ssemix.Update(bit,mixsse_upd_rate);
}

void BitplaneCoder::Encode(int32_t *abuf)
{
  pabuf=abuf;
  for (bpn=maxbpn;bpn>=0;bpn--)  {
    state=0;
    for (sample=0;sample<numsamples;sample++) {
      pestimate=PredictLaplace();
      GetSigState(sample);
      int bit=(pabuf[sample]>>bpn)&1;
      if (sigst[0]) { // coef is significant, refine
        rc.EncodeBitOne(PredictSSE(PredictRef()),bit);
        UpdateRef(bit);
        UpdateSSE(bit);
      } else { // coef is insignificant
        rc.EncodeBitOne(PredictSSE(PredictSig()),bit);
        UpdateSig(bit);
        UpdateSSE(bit);
        if (bit) msb[sample]=bpn;
      }
    }
  }
}

void BitplaneCoder::Decode(int32_t *buf)
{
  int bit;
  pabuf=buf;
  for (int i=0;i<numsamples;i++) buf[i]=0;
  for (bpn=maxbpn;bpn>=0;bpn--)  {
    state=0;
    for (sample=0;sample<numsamples;sample++) {
      pestimate=PredictLaplace();
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

