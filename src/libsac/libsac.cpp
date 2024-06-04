#include <algorithm>
#include "libsac.h"
#include "pred.h"
#include "../common/timer.h"
#include "../opt/dds.h"

FrameCoder::FrameCoder(int numchannels,int framesize,const coder_ctx &opt)
:pr_stages(numchannels, framesize),numchannels_(numchannels),framesize_(framesize),opt(opt)
{
  if (opt.profile==0) LoadProfileNormal(base_profile);
  else if (opt.profile==1) LoadProfileHigh(base_profile);

  profile_size_bytes_=base_profile.coefs.size()*4;

  framestats.resize(numchannels);
  samples.resize(numchannels);
  err0.resize(numchannels);
  err1.resize(numchannels);
  error.resize(numchannels);
  s2u_error.resize(numchannels);
  s2u_error_map.resize(numchannels);
  pred.resize(numchannels);
  for (int i=0;i<numchannels;i++) {
    samples[i].resize(framesize);
    err0[i].resize(framesize);
    err1[i].resize(framesize);
    error[i].resize(framesize);
    pred[i].resize(framesize);
    s2u_error[i].resize(framesize);
    s2u_error_map[i].resize(framesize);
  }
  encoded.resize(numchannels);
  enc_temp1.resize(numchannels);
  enc_temp2.resize(numchannels);
  numsamples_=0;
}

void FrameCoder::AnalyseChannel(int ch,int numsamples)
{
  if (numsamples) {
    int nbits=16;
    int minval=(1<<(nbits-1))-1,maxval=-(1<<(nbits-1));
    int64_t mean=0;
    int32_t *src=&(samples[ch][0]);

    uint32_t ordata=0,xordata=0,anddata=~0;

    for (int i=0;i<numsamples;i++) {
      if (src[i]<minval) minval=src[i];
      if (src[i]>maxval) maxval=src[i];

      xordata |= src[i] ^ -(src[i] & 1);
      anddata &= src[i];
      ordata |= src[i];

      mean+=src[i];
    }
    int tshift=0;
    if (!(ordata&1)) {
        while (!(ordata&1)) {
            tshift++;
            ordata>>=1;
        }
    } else if (anddata&1) {
        while (anddata&1) {
            tshift++;
            anddata>>=1;
        }
    } else if (!(xordata&2)) {
        while (!(xordata&2)) {
            tshift++;
            xordata>>=1;
        }
    }
    if (tshift) std::cout << "total shift: " << tshift << std::endl;

    mean=mean/static_cast<int64_t>(numsamples);
    std::cout << "mean: " << mean << '\n';
    if (mean) {
       minval-=mean;
       maxval-=mean;
       for (int i=0;i<numsamples;i++) src[i]-=mean;
    }
    framestats[ch].minval=minval;
    framestats[ch].maxval=maxval;
    framestats[ch].mean=mean;
    //std::cout << "mean: " << mean << ", min: " << minval << " ,max: " << maxval << std::endl;
  }
}


void SetParam(Predictor::tparam &param,const SacProfile &profile,bool optimize=false)
{
  param.nA=16;
  //param.nS0=round(profile.Get(9));
  //param.nS1=round(profile.Get(10));
  if (optimize) param.k=4;
  else param.k=1;
  param.lambda0=param.lambda1=profile.Get(0);
  param.ols_nu0=param.ols_nu1=profile.Get(1);
  param.mix_nu0=param.mix_nu1=1.0;
  param.bias_scale=32;
  param.bias_mu0 = param.bias_mu1 = 0.0015;
  //param.lpc_scale=0.0;

  if (profile.type==0) {
    //param.vn={256,32,4};
    param.vn0=param.vn1={256,32,4};
    param.vmu0={profile.Get(2),profile.Get(3),profile.Get(4)};
    param.vpowdecay0={profile.Get(5),profile.Get(6),profile.Get(7)};
    param.vmudecay0={profile.Get(8),profile.Get(9),profile.Get(10)};
    param.mu_mix0=profile.Get(11);
    param.mu_mix_beta0=profile.Get(12);

    param.lambda1=profile.Get(13);
    param.ols_nu1=profile.Get(14);
    param.vmu1={profile.Get(15),profile.Get(16),profile.Get(17)};
    param.vpowdecay1={profile.Get(18),profile.Get(19),profile.Get(20)};
    param.vmudecay1={profile.Get(21),profile.Get(22),profile.Get(23)};
    param.mu_mix1=profile.Get(24);
    param.mu_mix_beta1=profile.Get(25);

    param.nA=round(profile.Get(26));
    param.nB=round(profile.Get(27));
    param.nS0=round(profile.Get(28));
    param.nS1=round(profile.Get(29));

    param.beta_sum0=param.beta_sum1=profile.Get(30);
    param.beta_pow0=param.beta_pow1=profile.Get(31);
    param.beta_add0=param.beta_add1=profile.Get(32);

    param.bias_scale=profile.Get(33);

    param.mix_nu0=param.mix_nu1=profile.Get(34);

  } else if (profile.type==1) {
    param.vn0={(int)round(profile.Get(28)),(int)round(profile.Get(29)),(int)round(profile.Get(30)),(int)round(profile.Get(37))};
    param.vn1={(int)round(profile.Get(31)),(int)round(profile.Get(32)),(int)round(profile.Get(33)),(int)round(profile.Get(38))};
    //param.vn1={(int)round(profile.Get(28)),(int)round(profile.Get(29)),(int)round(profile.Get(30)),(int)round(profile.Get(37))};

    param.vmu0={profile.Get(2)/double(param.vn0[0]),profile.Get(3)/double(param.vn0[1]),profile.Get(4)/double(param.vn0[2]),profile.Get(5)/double(param.vn0[3])};
    param.vmudecay0={profile.Get(6),profile.Get(39),profile.Get(46),profile.Get(47)};
    param.vpowdecay0={profile.Get(7),profile.Get(8),profile.Get(50),profile.Get(51)};
    param.mu_mix0=profile.Get(10);
    param.mu_mix_beta0=profile.Get(11);

    param.lambda1=profile.Get(12);
    param.ols_nu1=profile.Get(13);
    param.vmu1={profile.Get(14)/double(param.vn1[0]),profile.Get(15)/double(param.vn1[1]),profile.Get(16)/double(param.vn1[2]),profile.Get(17)/double(param.vn1[3])};
    param.vmudecay1={profile.Get(18),profile.Get(40),profile.Get(48),profile.Get(49)};
    param.vpowdecay1={profile.Get(19),profile.Get(20),profile.Get(21),profile.Get(52)};
    param.mu_mix1=profile.Get(22);
    param.mu_mix_beta1=profile.Get(23);

    param.nA=round(profile.Get(24));
    param.nB=round(profile.Get(25));
    param.nS0=round(profile.Get(26));
    param.nS1=round(profile.Get(27));

    param.beta_sum0=profile.Get(34);
    param.beta_pow0=profile.Get(35);
    param.beta_add0=profile.Get(36);

    param.beta_sum1=profile.Get(34);
    param.beta_pow1=profile.Get(35);
    param.beta_add1=profile.Get(36);

    param.mix_nu0=profile.Get(41);
    param.mix_nu1=profile.Get(42);

    param.bias_mu0=profile.Get(43);
    param.bias_mu1=profile.Get(44);

    param.bias_scale=round(profile.Get(45));
  }
}

void FrameCoder::RemapError(int ch, int numsamples)
{
    int32_t *p=&(pred[ch][0]);
    int32_t *e=&(error[ch][0]);
    int32_t emax_map=0;

    for (int i=0;i<numsamples;i++) {
      int32_t map_e=framestats[ch].mymap.Map(p[i],e[i]);
      int32_t map_es=MathUtils::S2U(map_e);
      s2u_error_map[ch][i]=map_es;
      if (map_es>emax_map) emax_map=map_es;
    }
    framestats[ch].maxbpn_map=MathUtils::iLog2(emax_map);
}

/*double Cost_L1(const vec1D &src_e, int numsamples)
{
  double sum=0.0;
  if (numsamples) {
    for (int i=0;i<numsamples;i++)
    {
      sum+=std::fabs(src_e[i]);
    }
    sum /= double(numsamples);
  }
  return sum;
}*/

void FrameCoder::PredictStereoFrame(const SacProfile &profile,int ch0,int ch1,int from,int numsamples,bool optimize)
{
  const int32_t *src0=&(samples[ch0][from]);
  const int32_t *src1=&(samples[ch1][from]);

  int32_t *pred0=&(pred[ch0][0]);
  int32_t *dst0=&(error[ch0][0]);
  int32_t *pred1=&(pred[ch1][0]);
  int32_t *dst1=&(error[ch1][0]);

  Predictor::tparam param;
  SetParam(param,profile,optimize);
  Predictor pr(param);

  #if 0
  vec1D &stage0_ch0=pr_stages.pr_stage[ch0][0];
  vec1D &stage0_ch1=pr_stages.pr_stage[ch1][0];
  vec1D &stage1_ch0=pr_stages.pr_stage[ch0][1];
  vec1D &stage1_ch1=pr_stages.pr_stage[ch1][1];
  vec1D &stage2_ch0=pr_stages.pr_stage[ch0][2];
  vec1D &stage2_ch1=pr_stages.pr_stage[ch1][2];

  // stage0 ch0 OLS
  for (int i=0;i<numsamples;i++) {
    const double p=pr.Predict_stage0_ch0();
    const double val = src0[i];
    stage0_ch0[i]=val-p;
    pr.Update_stage0_ch0(val);
  }
  // stage0 ch1 stereo-OLS
  for (int i=0;i<numsamples;i++) {
    const double p=pr.Predict_stage0_ch1(src0,i,numsamples);
    const double val = src1[i];
    stage0_ch1[i]=val-p;
    pr.Update_stage0_ch1(val);
  }

  //stage1 ch0 NLMS
  for (int i=0;i<numsamples;i++) {
    const double p=pr.lms0.Predict();
    const double val = stage0_ch0[i];
    stage1_ch0[i]=val-p;
    pr.lms0.Update(val);
  }

  //stage1 ch1 NLMS
  for (int i=0;i<numsamples;i++) {
    const double p=pr.lms1.Predict();
    const double val = stage0_ch1[i];
    stage1_ch1[i]=val-p;
    pr.lms1.Update(val);
  }


  //if (!optimize) {
  //stage2 ch0 BC
  {
  for (int i=0;i<numsamples;i++) {
    const double val = src0[i];
    const double p_t = val - stage1_ch0[i];// calc p of OLS+NLMS
    const double p=pr.be0.Predict(p_t);
    stage2_ch0[i]=val-p;
    pr.be0.Update(val);
  }

  //stage2 ch1 BC
  for (int i=0;i<numsamples;i++) {
    const double val = src1[i];
    const double p_t = val - stage1_ch1[i];// calc p of OLS+NLMS
    const double p=pr.be1.Predict(p_t);
    stage2_ch1[i]=val-p;
  }
  }

  double *esrc0,*esrc1;
  if (optimize) {
    esrc0=&stage2_ch0[0];
    esrc1=&stage2_ch1[0];
  } else {
    esrc0=&stage2_ch0[0];
    esrc1=&stage2_ch1[0];
  }
  //store result
  for (int i=0;i<numsamples;i++) {
    const double val = src0[i];
    const double p = val - esrc0[i];
    const int32_t pi=clamp((int)std::round(p),framestats[0].minval,framestats[0].maxval);
    pred0[i] = pi;
    dst0[i]=val-pi;
  }
  for (int i=0;i<numsamples;i++) {
    const double val = src1[i];
    const double p = val - esrc1[i];
    const int32_t pi=clamp((int)std::round(p),framestats[1].minval,framestats[1].maxval);
    pred1[i] = pi;
    dst1[i]=val-pi;
  }

  #if 0
    std::cout << Cost_L1(stage0_ch0,numsamples) << '\n';
    std::cout << Cost_L1(stage1_ch0,numsamples) << '\n';
    std::cout << Cost_L1(stage2_ch0,numsamples) << '\n';
  #endif

  #else
  for (int i=0;i<numsamples;i++) {
    double p=pr.PredictMaster();

    int32_t p0=clamp((int)std::round(p),framestats[0].minval,framestats[0].maxval);
    pred0[i]=p0;
    dst0[i]=src0[i]-p0;

    pr.UpdateMaster(src0[i]);
  }
  for (int i=0;i<numsamples;i++) {
    double p=pr.PredictSlave(src0,i,numsamples);//lm1.Predict();
    int32_t p1=clamp((int)std::round(p),framestats[1].minval,framestats[1].maxval);

    pred1[i]=p1;
    dst1[i]=src1[i]-p1;

    pr.UpdateSlave(src1[i]);
  }
  #endif

  int32_t emax0=0,emax1=0;
  for (int i=0;i<numsamples;i++) {
    int32_t e0=MathUtils::S2U(dst0[i]);
    if (e0>emax0) emax0=e0;
    s2u_error[ch0][i]=e0;

    int32_t e1=MathUtils::S2U(dst1[i]);
    if (e1>emax1) emax1=e1;
    s2u_error[ch1][i]=e1;
  }
  framestats[0].maxbpn=MathUtils::iLog2(emax0);
  framestats[1].maxbpn=MathUtils::iLog2(emax1);
}

void FrameCoder::UnpredictStereoFrame(const SacProfile &profile,int ch0,int ch1,int numsamples)
{
  const int32_t *src0=&(error[ch0][0]);
  const int32_t *src1=&(error[ch1][0]);
  int32_t *dst0=&(samples[ch0][0]);
  int32_t *dst1=&(samples[ch1][0]);

  Predictor::tparam param;
  SetParam(param,profile,false);
  Predictor pr(param);
  for (int i=0;i<numsamples;i++) {
    double p=pr.PredictMaster();
    int32_t p0=clamp((int)round(p),framestats[0].minval,framestats[0].maxval);

    if (framestats[0].enc_mapped) dst0[i]=p0+framestats[0].mymap.Unmap(p0,src0[i]);
    else dst0[i]=p0+src0[i];

    pr.UpdateMaster(dst0[i]);
  }
  for (int i=0;i<numsamples;i++) {
    double p=pr.PredictSlave(dst0,i,numsamples);
    int32_t p1=clamp((int)round(p),framestats[1].minval,framestats[1].maxval);

    if (framestats[1].enc_mapped) dst1[i]=p1+framestats[1].mymap.Unmap(p1,src1[i]);
    else dst1[i]=p1+src1[i];

    pr.UpdateSlave(dst1[i]);
  }
}

int FrameCoder::EncodeMonoFrame_Normal(int ch,int numsamples,BufIO &buf)
{
  buf.Reset();
  RangeCoderSH rc(buf);
  rc.Init();
  BitplaneCoder bc(rc,framestats[ch].maxbpn,numsamples);
  int32_t *psrc=&(s2u_error[ch][0]);
  bc.Encode(psrc);
  rc.Stop();
  return buf.GetBufPos();
}

int FrameCoder::EncodeMonoFrame_Mapped(int ch,int numsamples,BufIO &buf)
{
  buf.Reset();

  RangeCoderSH rc(buf);
  rc.Init();

  BitplaneCoder bc(rc,framestats[ch].maxbpn_map,numsamples);

  MapEncoder me(rc,framestats[ch].mymap.usedl,framestats[ch].mymap.usedh);
  me.Encode();
  //std::cout << "mapsize: " << buf.GetBufPos() << " Bytes\n";

  bc.Encode(&(s2u_error_map[ch][0]));
  rc.Stop();
  return buf.GetBufPos();
}

void FrameCoder::EncodeMonoFrame(int ch,int numsamples)
{
  if (opt.sparse_pcm==0) {
    EncodeMonoFrame_Normal(ch,numsamples,enc_temp1[ch]);
    framestats[ch].enc_mapped=false;
    encoded[ch]=enc_temp1[ch];
  } else {
    RemapError(0,numsamples);
    RemapError(1,numsamples);
    int size_normal=EncodeMonoFrame_Normal(ch,numsamples,enc_temp1[ch]);
    int size_mapped=EncodeMonoFrame_Mapped(ch,numsamples,enc_temp2[ch]);
    if (opt.verbose_level>0) {
      std::cout << "block " << size_normal << ",map " << size_mapped << '\n';
    }
    if (size_normal<size_mapped)
    {
      //std::cout << "block: normal\n";
      framestats[ch].enc_mapped=false;
      encoded[ch]=enc_temp1[ch];
    } else {
      framestats[ch].enc_mapped=true;
      encoded[ch]=enc_temp2[ch];
    }
  }
}

void FrameCoder::DecodeMonoFrame(int ch,int numsamples)
{
  int32_t *dst=&(error[ch][0]);
  BufIO &buf=encoded[ch];
  buf.Reset();

  RangeCoderSH rc(buf,1);
  rc.Init();
  if (framestats[ch].enc_mapped) {
    framestats[ch].mymap.Reset();
    MapEncoder me(rc,framestats[ch].mymap.usedl,framestats[ch].mymap.usedh);
    me.Decode();
    //std::cout << buf.GetBufPos() << std::endl;
  }

  BitplaneCoder bc(rc,framestats[ch].maxbpn,numsamples);
  bc.Decode(dst);
  rc.Stop();
}

double FrameCoder::GetCost(SacProfile &profile,CostFunction *func,int start_sample,int samples_to_optimize)
{
  const int32_t *err0=&(error[0][0]);
  const int32_t *err1=&(error[1][0]);

  PredictStereoFrame(profile,0,1,start_sample,samples_to_optimize,true);
  double c1=func->Calc(err0,samples_to_optimize);
  double c2=func->Calc(err1,samples_to_optimize);
  return c1+c2;
}

void PrintProfile(int iprofile, SacProfile &profile)
{
    Predictor::tparam param;
    SetParam(param,profile);

    if (iprofile==0) {
      std::cout << "lpc ";
      std::cout << "nA " << std::round(profile.Get(26)) << ' ' << "nB " << std::round(profile.Get(27)) << ' ';
      std::cout << "nS0 " << std::round(profile.Get(28)) << ' ' << "nS1 " << std::round(profile.Get(29)) << '\n';
    } else if (iprofile==1) {
      std::cout << '\n';
      std::cout << "lpc ";
      std::cout << "nA " << std::round(profile.Get(24)) << ' ' << "nB " << std::round(profile.Get(25)) << ' ';
      std::cout << "nS0 " << std::round(profile.Get(26)) << ' ' << "nS1 " << std::round(profile.Get(27)) << '\n';
      std::cout << "lms0 ";
      for (int i=28;i<=30;i++) std::cout << std::round(profile.Get(i)) << ' ';
      std::cout << std::round(profile.Get(37));
      std::cout << '\n';
      std::cout << "lms1 ";
      for (int i=31;i<=33;i++) std::cout << std::round(profile.Get(i)) << ' ';
      std::cout << std::round(profile.Get(38));
      std::cout << '\n';
      std::cout << "mu_decay ";
      for (auto &x : param.vmudecay0)
        std::cout << x << ' ';
      std::cout << '\n';
      std::cout << "pow_decay ";
      for (auto &x : param.vpowdecay0)
        std::cout << x << ' ';
      std::cout << '\n';
    }
    std::cout << "mu-nu " << param.mix_nu0 << ", " << param.mix_nu1 << "\n";
    std::cout << "bias mu " << param.bias_mu0 << ", " << param.bias_mu1 << ", scale " << param.bias_scale << "\n";
    std::cout << "lpc nu " << param.ols_nu0 << ' ' << param.ols_nu1 << '\n';
    std::cout << "lpc cov0 " << param.beta_sum0 << ' ' << param.beta_pow0 << ' ' << param.beta_add0 << "\n";
}

void FrameCoder::Optimize(SacProfile &profile,const std::vector<int>&params_to_optimize)
{
  //std::cout << params_to_optimize.size() << " parameters to optimize\n";
  int norm_samples=std::min((int)std::ceil(framesize_*opt.optimize_fraction), framesize_);
  int samples_to_optimize=std::min((int)std::ceil(numsamples_*opt.optimize_fraction), framesize_);
  if (samples_to_optimize<norm_samples) {
    samples_to_optimize=std::min(numsamples_,norm_samples);
  }

  const int start_pos=(numsamples_-samples_to_optimize)/2;

  CostFunction *CostFunc=nullptr;
  switch (opt.optimize_cost)  {
    case opt.SearchCost::L1:CostFunc=new CostMeanRMS();break;
    case opt.SearchCost::Golomb:CostFunc=new CostGolomb();break;
    case opt.SearchCost::Entropy:CostFunc=new CostEntropyO0();break;
    case opt.SearchCost::Bitplane:CostFunc=new CostBitplane();break;
    default:std::cerr << "  error: unknown FramerCoder::CostFunction\n";return;
  }

  if (opt.optimize_search==opt.SearchMethod::DDS) {
    const int ndim=params_to_optimize.size();
    vec1D xstart(ndim); // starting vector
    Opt::box_const pb(ndim); // set constraints
    for (int i=0;i<ndim;i++) {
      pb[i].xmin=profile.coefs[params_to_optimize[i]].vmin;
      pb[i].xmax=profile.coefs[params_to_optimize[i]].vmax;
      xstart[i]=profile.coefs[params_to_optimize[i]].vdef;
      //std::cout << pb[i].xmin << ' ' << pb[i].xmax << ' ' << xstart[i] << '\n';
    }

    auto cost_func=[&](const vec1D &x) {
      for (int i=0;i<ndim;i++) profile.coefs[params_to_optimize[i]].vdef=x[i];
      return GetCost(profile,CostFunc,start_pos,samples_to_optimize);
    };


    //Opt::opt_ret ret;

    DDS dds(pb);

    if (opt.verbose_level>0) std::cout << "\nDDS " << opt.optimize_maxnfunc << "= ";
    auto ret = dds.run(cost_func,xstart,opt.optimize_maxnfunc,opt.dds_search_radius);
    if (opt.verbose_level>0) std::cout << ret.first << "\n";

    const vec1D x_p = ret.second;
    for (int i=0;i<ndim;i++) profile.coefs[params_to_optimize[i]].vdef=x_p[i]; // save optimal vector
  }

  if (opt.verbose_level>0) {
    PrintProfile(opt.profile, profile);
  }
  delete CostFunc;
}


void Codec::SetOptimizeParam(FrameCoder::coder_ctx &opt)
{
  opt.optimize_fraction=0.25;
  opt.dds_search_radius=0.15;
  opt.optimize_cost=opt.SearchCost::Entropy;
  opt.optimize_search=opt.SearchMethod::DDS;

  if (opt.optimize_mode==0) {
    opt.optimize_maxnfunc=150;
    opt.optimize_fraction=0.20;
  } else if (opt.optimize_mode==1) {
    opt.optimize_maxnfunc=500;
    opt.optimize_fraction=0.20;
  } else if (opt.optimize_mode==2) {
    opt.optimize_maxnfunc=1000;
    opt.optimize_fraction=0.50;
  } else if (opt.optimize_mode==3) {
    opt.optimize_maxnfunc=1500;
    opt.optimize_fraction=0.70;
  } else if (opt.optimize_mode==4) { // insane
    opt.optimize_cost=opt.SearchCost::Bitplane;
    opt.optimize_maxnfunc=2000;
    opt.optimize_fraction=1.0;
  } else if (opt.optimize_mode==5) { // bias-param
    opt.optimize_maxnfunc=100;
  }
}

double FrameCoder::AnalyseStereoChannel(int ch0, int ch1, int numsamples)
{
  int32_t *src0=&(samples[ch0][0]);
  int32_t *src1=&(samples[ch1][0]);
  int64_t sum0=0,sum1=0,sum_m=0,sum_s=0;
  for (int i=0;i<numsamples;i++) {
    sum0+=fabs(src0[i]);
    sum1+=fabs(src1[i]);
    int32_t m=(src0[i]+src1[i]) / 2;
    int32_t s=(src0[i]-src1[i]);

    sum_m+=fabs(m);
    sum_s+=fabs(s);
  }
  int64_t c0 = sum0+sum1;
  int64_t c1 = sum_m+sum_s;
  return double(c0) / double(c1);
}

void FrameCoder::ApplyMs(int ch0, int ch1, int numsamples)
{
  int32_t *src0=&(samples[ch0][0]);
  int32_t *src1=&(samples[ch1][0]);
  for (int i=0;i<numsamples;i++) {
    int32_t m=(src0[i]+src1[i]) / 2;
    int32_t s=(src0[i]-src1[i]);
    src0[i]=m;
    src1[i]=s;
  }
}

void FrameCoder::AnalyseMonoChannel(int ch, int numsamples)
{
  int32_t *src=&(samples[ch][0]);
  int64_t sum=0;

  if (numsamples) {
    for (int i=0;i<numsamples;i++) {
        sum += src[i];
    }
    framestats[ch].mean = (sum) / numsamples;

    int32_t minval = std::numeric_limits<int32_t>::max();
    int32_t maxval = std::numeric_limits<int32_t>::min();
    for (int i=0;i<numsamples;i++) {
      const int32_t val=src[i];
      if (val>maxval) maxval=val;
      else if (val<minval) minval=val;
    }
    framestats[ch].minval = minval;
    framestats[ch].maxval = maxval;
    //std::cout << "stats: " << numsamples << '\n';
    //std::cout << "mean:" << framestats[ch].mean << " min:" << framestats[ch].minval << " max:" << framestats[ch].maxval << "\n";
  }
}




void FrameCoder::Predict()
{
  for (int ch=0;ch<numchannels_;ch++)
  {
    //AnalyseChannel(ch,numsamples_);
    AnalyseMonoChannel(ch,numsamples_);
    framestats[ch].mymap.Reset();
    framestats[ch].mymap.Analyse(&(samples[ch][0]),numsamples_);
  }

  if (opt.stereo_ms==1 && numchannels_==2) {
     double r=AnalyseStereoChannel(0,1,numsamples_);
     if (opt.verbose_level>0) {
       std::cout << "ms-r " << r << '\n';
     }
     if (r > 1.5) {
       ApplyMs(0,1,numsamples_);
     }
  }

  if (opt.optimize)
  {
    // reset profile params
    // otherwise: starting point for optimization is the best point from the last frame
    if (opt.reset_profile) {
      if (opt.profile==0) LoadProfileNormal(base_profile);
      else if (opt.profile==1) LoadProfileHigh(base_profile);
    }
    if (opt.profile==0) {
       std::vector<int>lparam(base_profile.coefs.size());
       std::iota(std::begin(lparam),std::end(lparam),0);
       Optimize(base_profile,lparam);
    } else if (opt.profile==1) {
       std::vector<int>lparam(base_profile.coefs.size());
       std::iota(std::begin(lparam),std::end(lparam),0);
       Optimize(base_profile,lparam);
     }
  }
  if (numchannels_==2) {
    PredictStereoFrame(base_profile,0,1,0,numsamples_,false);
  }
}

void FrameCoder::Unpredict()
{
  if (numchannels_==2) UnpredictStereoFrame(base_profile,0,1,numsamples_);
  /*for (int ch=0;ch<numchannels_;ch++) {
    if (framestats[ch].mean != 0) {
      int32_t *dst=&(samples[ch][0]);
      for (int i=0;i<numsamples_;i++) dst[i] += framestats[ch].mean;
    }
  }*/
  //else for (int ch=0;ch<numchannels;ch++) UnpredictMonoFrame(ch,numsamples);
}

void FrameCoder::Encode()
{
  for (int ch=0;ch<numchannels_;ch++) EncodeMonoFrame(ch,numsamples_);
}

void FrameCoder::Decode()
{
  for (int ch=0;ch<numchannels_;ch++) DecodeMonoFrame(ch,numsamples_);
}

void FrameCoder::EncodeProfile(const SacProfile &profile,std::vector <uint8_t>&buf)
{
  //assert(sizeof(float)==4);
  //std::cout << "number of coefs: " << profile.coefs.size() << " (" << profile_size_bytes_ << ")" << std::endl;

  uint32_t ix;
  for (int i=0;i<(int)profile.coefs.size();i++) {
     memcpy(&ix,&profile.coefs[i].vdef,4);
     //ix=*((uint32_t*)&profile.coefs[i].vdef);
     BitUtils::put32LH(&buf[4*i],ix);
  }
}

void FrameCoder::DecodeProfile(SacProfile &profile,const std::vector <uint8_t>&buf)
{
  uint32_t ix;
  for (int i=0;i<(int)profile.coefs.size();i++) {
     ix=BitUtils::get32LH(&buf[4*i]);
     memcpy(&profile.coefs[i].vdef,&ix,4);
     //profile.coefs[i].vdef=*((float*)&ix);
  }
}

int FrameCoder::WriteBlockHeader(std::fstream &file, std::vector<SacProfile::FrameStats> &framestats,int ch)
{
  uint8_t buf[12];
  BitUtils::put32LH(buf,framestats[ch].blocksize);
  BitUtils::put16LH(buf+4,framestats[ch].mean);
  BitUtils::put16LH(buf+6,framestats[ch].minval);
  BitUtils::put16LH(buf+8,framestats[ch].maxval);
  uint16_t flag=0;
  if (framestats[ch].enc_mapped) {
     flag|=(1<<9);
     flag|=framestats[ch].maxbpn_map;
  } else {
     flag|=framestats[ch].maxbpn;
  }
  BitUtils::put16LH(buf+10,flag);
  file.write(reinterpret_cast<char*>(buf),12);
  return 12;
}

int FrameCoder::ReadBlockHeader(std::fstream &file, std::vector<SacProfile::FrameStats> &framestats,int ch)
{
  uint8_t buf[12];
  file.read(reinterpret_cast<char*>(buf),12);

  framestats[ch].blocksize=BitUtils::get32LH(buf);
  framestats[ch].mean=static_cast<int16_t>(BitUtils::get16LH(buf+4));
  framestats[ch].minval=static_cast<int16_t>(BitUtils::get16LH(buf+6));
  framestats[ch].maxval=static_cast<int16_t>(BitUtils::get16LH(buf+8));
  uint16_t flag=BitUtils::get16LH(buf+10);
  if (flag>>9) framestats[ch].enc_mapped=true;
  else framestats[ch].enc_mapped=false;
  framestats[ch].maxbpn=flag&0xff;
  return 12;
}

void FrameCoder::WriteEncoded(AudioFile &fout)
{
  uint8_t buf[12];
  BitUtils::put32LH(buf,numsamples_);
  fout.file.write(reinterpret_cast<char*>(buf),4);
  std::vector <uint8_t>profile_buf(profile_size_bytes_);
  EncodeProfile(base_profile,profile_buf);
  fout.file.write(reinterpret_cast<char*>(&profile_buf[0]),profile_size_bytes_);
  for (int ch=0;ch<numchannels_;ch++) {
    framestats[ch].blocksize = encoded[ch].GetBufPos();
    WriteBlockHeader(fout.file, framestats, ch);
    fout.WriteData(encoded[ch].GetBuf(),framestats[ch].blocksize);
  }
}

void FrameCoder::ReadEncoded(AudioFile &fin)
{
  uint8_t buf[8];
  fin.file.read(reinterpret_cast<char*>(buf),4);
  numsamples_=BitUtils::get32LH(buf);
  std::vector <uint8_t>profile_buf(profile_size_bytes_);
  fin.file.read(reinterpret_cast<char*>(&profile_buf[0]),profile_size_bytes_);
  DecodeProfile(base_profile,profile_buf);

  for (int ch=0;ch<numchannels_;ch++) {
    ReadBlockHeader(fin.file, framestats, ch);
    fin.ReadData(encoded[ch].GetBuf(),framestats[ch].blocksize);
  }
}

void Codec::PrintProgress(int samplesprocessed,int totalsamples)
{
  double r=samplesprocessed*100.0/(double)totalsamples;
  std::cout << "  " << samplesprocessed << "/" << totalsamples << ":" << std::setw(6) << miscUtils::ConvertFixed(r,1) << "%\r";
}

void Codec::ScanFrames(Sac &mySac)
{
  std::vector<SacProfile::FrameStats> framestats(mySac.getNumChannels());
  std::streampos fsize=mySac.getFileSize();

  SacProfile profile_tmp; //create dummy profile
  switch (mySac.GetProfile())
  {
    case 0:LoadProfileNormal(profile_tmp);break;
    case 1:LoadProfileHigh(profile_tmp);break;
  }
  const int size_profile_bytes=profile_tmp.coefs.size()*4;

  int frame_num=1;
  while (mySac.file.tellg()<fsize) {
    uint8_t buf[12];
    mySac.file.read(reinterpret_cast<char*>(buf),4);
    int numsamples=BitUtils::get32LH(buf);
    std::cout << "Frame " << frame_num << ": " << numsamples << " samples "<< std::endl;

    mySac.file.seekg(size_profile_bytes,std::ios_base::cur); // skip profile coefs

    for (int ch=0;ch<mySac.getNumChannels();ch++) {
      FrameCoder::ReadBlockHeader(mySac.file, framestats, ch);
      std::cout << "  Channel " << ch << ": " << framestats[ch].blocksize << " bytes\n";
      std::cout << "    Bpn: " << framestats[ch].maxbpn << ", sparse_pcm: " << (framestats[ch].enc_mapped) << std::endl;
      std::cout << "    mean: " << framestats[ch].mean << ", min: " << framestats[ch].minval << ", max: " << framestats[ch].maxval << std::endl;
      mySac.file.seekg(framestats[ch].blocksize, std::ios_base::cur);
    }
    frame_num++;
  }
  std::cout << "total frames: " << (frame_num-1) << '\n';
}

void Codec::EncodeFile(Wav &myWav,Sac &mySac,FrameCoder::coder_ctx &opt)
{
  const int numchannels=myWav.getNumChannels();
  framesize=8*myWav.getSampleRate();

  SetOptimizeParam(opt);
  FrameCoder myFrame(numchannels,framesize,opt);

  mySac.SetProfile(opt.profile);
  mySac.WriteHeader(myWav);
  std::streampos hdrpos = mySac.file.tellg();
  mySac.WriteMD5(myWav.md5ctx.digest);
  myWav.InitFileBuf(framesize);

  Timer gtimer,ltimer;
  double time_prd,time_enc;
  time_prd=time_enc=0.;

  gtimer.start();
  int samplescoded=0;
  int samplestocode=myWav.getNumSamples();
  while (samplestocode>0) {
    int samplesread=myWav.ReadSamples(myFrame.samples,framesize);

    myFrame.SetNumSamples(samplesread);

    ltimer.start();myFrame.Predict();ltimer.stop();time_prd+=ltimer.elapsedS();
    ltimer.start();myFrame.Encode();ltimer.stop();time_enc+=ltimer.elapsedS();
    myFrame.WriteEncoded(mySac);

    samplescoded+=samplesread;
    PrintProgress(samplescoded,myWav.getNumSamples());
    samplestocode-=samplesread;
  }
  MD5::Finalize(&myWav.md5ctx);
  gtimer.stop();
  double time_total=gtimer.elapsedS();
  if (time_total>0.)   {
     double rprd=time_prd*100./time_total;
     double renc=time_enc*100./time_total;
     std::cout << "  Timing: pred " << miscUtils::ConvertFixed(rprd,2) << "%, ";
     std::cout << "enc " << miscUtils::ConvertFixed(renc,2) << "%, ";
     std::cout << "misc " << miscUtils::ConvertFixed(100.-rprd-renc,2) << "%" << std::endl;
  }
  std::cout << "\n  Audio MD5: ";
  for (auto x : myWav.md5ctx.digest) std::cout << std::hex << (int)x;
  std::cout << std::dec << '\n';

  std::streampos eofpos = mySac.file.tellg();
  mySac.file.seekg(hdrpos);
  mySac.WriteMD5(myWav.md5ctx.digest);
  mySac.file.seekg(eofpos);
}

void Codec::DecodeFile(Sac &mySac,Wav &myWav)
{
  const int numchannels=mySac.getNumChannels();
  framesize=8*mySac.getSampleRate();

  int samplestodecode=mySac.getNumSamples();
  int samplesdecoded=0;

  myWav.InitFileBuf(framesize);
  mySac.UnpackMetaData(myWav);
  myWav.WriteHeader();
  FrameCoder::coder_ctx opt;
  opt.profile=mySac.GetProfile();
  FrameCoder myFrame(numchannels,framesize,opt);
  while (samplestodecode>0) {
    myFrame.ReadEncoded(mySac);
    myFrame.Decode();
    myFrame.Unpredict();
    myWav.WriteSamples(myFrame.samples,myFrame.GetNumSamples());

    samplesdecoded+=myFrame.GetNumSamples();
    PrintProgress(samplesdecoded,myWav.getNumSamples());
    samplestodecode-=myFrame.GetNumSamples();
  }
  myWav.WriteHeader();
}
