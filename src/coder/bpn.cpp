#include "../common/utils.h"
#include "bpn.h"

BitplaneCoder::BitplaneCoder(int maxbpn,int numsamples)
:csig0(1<<20),csig1(1<<20),csig2(1<<20),
cref0(1<<20),cref1(1<<20),cref2(1<<20),
lmixref(256,LogMixer(5)),lmixsig(256,LogMixer(4)),
ssemix(2),
sm(4,ScaleModel()),
sigstats(SIG_RADIUS,SIG_MAXAGE),
msb(numsamples),
maxbpn(maxbpn),numsamples(numsamples)
{
  bpn=0;
  nrun=0;
  for (int i=0;i<32;i++)
    bmask[i]=~((uint32_t{1}<<i)-1);
  //set initial mixer weights
  for (int i=0;i<256;i++) {
    lmixref[i].SetConfidence(0.01);
    //set laplace weight higher
    lmixref[i].w[0] = WSCALE/2;
    lmixsig[i].SetConfidence(0.15);
  }
  //create quantization table for activity context
  int num_qcontext=sigstats.get_max_acc()+1;
  q_logbin.resize(num_qcontext);
  for (int i=0;i<num_qcontext;i++) {
    q_logbin[i]=LogBin(i,1.0,1.5,32);
  }
}

// x is a fixed-point logit: real logit = x / LSCALE.
int BitplaneCoder::LaplaceLogitBin(int x, int nbits)
{
  assert(nbits >= 0 && nbits <= 16);
  const int lo = BM::DMIN;
  const int hi = 0; //p_laplace is always <= 0.5
  const int nbins=1<<nbits;
  x = std::clamp(x, lo, hi);
  return std::min(nbins-1,((x-lo)*nbins)/(hi - lo));
}

int BitplaneCoder::LogBin(unsigned x,double k,double g,int nbins)
{
  if (x == 0) return 0;
  const double y = std::log1p(static_cast<double>(x)*k)/std::log(g);
  return static_cast<int>(std::clamp(std::ceil(y),1.0,double(nbins-1)));
}

//fast LogBin with k=1.0,g=2.0
int BitplaneCoder::LogBin(unsigned x,int nbins)
{
  return std::min(BitUtils::count_bits32(x),nbins-1);
}

static int RelativeMagnitudeClass(double r)
{
  return r < 0.25 ? 0:
         r < 0.5  ? 1:
         r < 1.0  ? 2:
         r < 2.0  ? 3:
         r < 4.0  ? 4:
         r < 8.0  ? 5:
         r < 16.0 ? 6:7;
}

void BitplaneCoder::GetState(double mean,int sample,State &s)
{
  int *x=s.sigst;
  x[0]=msb[sample];
  const int maxr=8;
  const int maxr_age=4;
  const int maxr_agecap=8;
  s.age_sum=0;
  for (int r=1;r<=maxr;r++)
  {
    int &xl=x[2*r-1];
    int &xr=x[2*r];
    xl = sample>=r?msb[sample-r]:0;
    xr = sample<numsamples-r?msb[sample+r]:0;

    if (r<=maxr_age) {
      s.age_sum += xl?std::min(xl-bpn,maxr_agecap):0;
      s.age_sum += xr?std::min(xr-bpn,maxr_agecap):0;
    }
  }

  const double rho=0.5;
  const int plane= uint32_t{1}<<bpn;
  const int value = static_cast<uint32_t>(pabuf[sample]);
  const int prefix = (value >> (bpn+1)) << (bpn+1);
  const double rel_ratio = static_cast<double>(prefix) /(mean+rho*static_cast<double>(plane));
  s.mclass = RelativeMagnitudeClass(rel_ratio);
}

int BitplaneCoder::PredictRef()
{
  int val=pabuf[sample];

  int lval=sample>0?pabuf[sample-1]:0;
  int lval2=sample>1?pabuf[sample-2]:0;
  int nval=sample<(numsamples-1)?pabuf[sample+1]:0;
  int nval2=sample<(numsamples-2)?pabuf[sample+2]:0;

  int shift0=bpn;
  int shift1=bpn+1;

  int b0=val>>shift1;
  int b1=lval>>shift0;
  int b2=nval>>shift1;
  int b3=lval2>>shift0;
  int b4=nval2>>shift1;

  int x0=b0<<1;
  int x1=b1;
  int x2=b2<<1;
  int x3=b3;
  int x4=b4<<1;
  int d0=4*x0>(x1+x2+x3+x4);

  int c0=x0<x1?1:0;
  int c1=x0<x2?1:0;
  int c2=x0<x3?1:0;
  int c3=x0<x4?1:0;

  int ctx0=s.mclass|(std::min(s.age_sum,31)<<3);//8-bit
  int ctx1=(b0&31)|((b1&7)<<5)|((b3&1)<<8); //9-bits
  int ctx2=c0|(c1<<1)|(c2<<2)|(c3<<3)|(d0<<4); //5-bits

  pc1=&cref0[ctx0];
  pc2=&cref1[ctx1];
  pc3=&cref2[ctx2];

  //8-bit mixer context
  int ref_mixctx=(std::min(msb[sample]-bpn,15)<<4)|(s.mclass<<1)|d0;
  plmix=&lmixref[ref_mixctx];

  const int input[] = {s.pl_prob,pl.p1,pc1->p1,pc2->p1,pc3->p1};
  return plmix->Predict(input);
}

void BitplaneCoder::UpdateRef(int bit)
{
  pl.update(bit,cnt_upd_rate_p);
  pc1->update(bit,cnt_upd_rate_ref);
  pc2->update(bit,cnt_upd_rate_ref);
  pc3->update(bit,cnt_upd_rate_ref);
  plmix->Update(bit,mix_upd_rate_ref);
}

// for testing: count number of significant samples in neighborhood
std::tuple<int,int,int> BitplaneCoder::CountSig(int n)
{
  int n1=0,n2=0,n3=0;

  auto update_scores=[&](int x,int i)
  {
    if (x) {
      const int d=std::min(x-bpn,2);
      n1+=1;
      if (d) {
        n2++;
        n3+=i*d;
      }
    }
  };
  for (int i=1;i<=n;i++) {
    if (sample-i>=0)
      update_scores(msb[sample-i],i);

    if (sample+i<numsamples)
      update_scores(msb[sample+i],i);
  }
  return {n1,n2,n3};
}

int BitplaneCoder::PredictSig()
{
  int ctx1=0;
  for (int i=0;i<16;i++)
    if (s.sigst[i+1]) ctx1|=1<<i;

  const SigStats::St c = sigstats.Counts();
  assert(c.n4 < q_logbin.size());
  const int ctx2=(c.n2<<5)|q_logbin[c.n4];

  pc1=&csig0[ctx1];
  pc2=&csig1[ctx2];

  //8-bit
  int sig_mixctx=(nrun<<4)|((c.n1>=3?3:c.n1)<<2)|((c.n2>0)<<1)|(c.n1-c.n2>1);
  plmix=&lmixsig[sig_mixctx];
  const int input[]={s.pl_scale,pl.p1,pc1->p1,pc2->p1};
  return plmix->Predict(input);
}

void BitplaneCoder::UpdateSig(int bit)
{
  pl.update(bit,cnt_upd_rate_p);
  pc1->update(bit,cnt_upd_rate_sig);
  pc2->update(bit,cnt_upd_rate_sig);
  plmix->Update(bit,mix_upd_rate_sig);
}

int BitplaneCoder::PredictSSE(int p1)
{
  int dist=s.sigst[0]?std::min(msb[sample]-bpn,7):0;
  int ctx1=(LaplaceLogitBin(s.pl_logit,5)<<3)+(dist); //8-bit
  int ctx2=(s.sigst[0]?1:0)+((s.sigst[1]?1:0)<<1)+((s.sigst[2]?1:0)<<2)+(std::min(s.age_sum,31)<<3); //8-bit
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

template <bool IsDecode,class Coder>
void BitplaneCoder::Run(Coder &&code,int32_t *buf)
{
  pabuf=buf;
  ExpMean em(0.95,numsamples);
  LaplaceModel lm;
  if constexpr(IsDecode)
    std::fill_n(buf,numsamples,0);

  for (bpn=maxbpn;bpn>=0;bpn--)  {
    const uint32_t plane=uint32_t{1}<<bpn;
    em.Begin(pabuf,bpn);
    for (auto& model:sm)
      model.Begin();
    sigstats.Begin(msb,bpn);
    nrun=0;

    for (sample=0;sample<numsamples;sample++) {
      const double mean=em.Mean();
      s.pl_logit=lm.PredictLogit(mean,plane);
      s.pl_prob=myDomain.Inv(s.pl_logit);
      if (sample==0) pl.p1 = s.pl_prob;

      GetState(mean,sample,s);

      int bit=0;
      if constexpr(!IsDecode) // input bit
        bit=(pabuf[sample]>>bpn)&1;

      if (s.sigst[0]) { // coef is significant, refine
        bit = code(PredictSSE(PredictRef()),bit);
        UpdateRef(bit);
        UpdateSSE(bit);
        nrun=0;
      } else { // coef is insignificant
        const int sctx=nrun>8?2:nrun>0?1:0;
        s.pl_scale = sm[sctx].Predict(mean,plane);
        bit = code(PredictSSE(PredictSig()),bit);
        UpdateSig(bit);
        UpdateSSE(bit);
        sm[sctx].Update(bit);
        if (bit) msb[sample]=bpn;
        if (nrun<15) nrun++;
      }
      if constexpr(IsDecode)
        buf[sample] |= bit << bpn;

      if (sample<numsamples-1) {
        em.Advance();
        sigstats.Advance();
      }
    }
  }
  if constexpr(IsDecode)
    for (int i=0;i<numsamples;i++)
      buf[i]=MathUtils::U2S(buf[i]);
}

void BitplaneCoder::Encode(EncodeP1 encode_p1,int32_t *buf)
{
  Run<false>(
    [&](uint32_t p1, int bit) -> int {
      encode_p1(p1, bit);
      return bit;
    },buf);
}

void BitplaneCoder::Decode(DecodeP1 decode_p1,int32_t *buf)
{
  Run<true>(
    [&](uint32_t p1, int) -> int {
      return decode_p1(p1);
    },buf);
}

