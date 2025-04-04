#include <algorithm>
#include <cmath>
#include <thread>
#include <future>

#include "libsac.h"
#include "pred.h"
#include "sparse.h"
#include "../common/timer.h"
#include "../opt/dds.h"
#include "../opt/de.h"

FrameCoder::FrameCoder(int numchannels,int framesize,const coder_ctx &opt)
:numchannels_(numchannels),framesize_(framesize),opt(opt)
{
  profile_size_bytes_=base_profile.LoadBaseProfile()*4;

  framestats.resize(numchannels);
  samples.resize(numchannels);
  error.resize(numchannels);
  s2u_error.resize(numchannels);
  s2u_error_map.resize(numchannels);
  pred.resize(numchannels);
  for (int i=0;i<numchannels;i++) {
    samples[i].resize(framesize);
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



void FrameCoder::SetParam(Predictor::tparam &param,const SacProfile &profile,bool optimize)
{
  if (optimize) param.k=opt.ocfg.optk;
  else param.k=1;

  param.lambda0=param.lambda1=profile.Get(0);
  param.ols_nu0=param.ols_nu1=profile.Get(1);

  param.vn0={(int)round(profile.Get(28)),(int)round(profile.Get(29)),(int)round(profile.Get(30)),(int)round(profile.Get(37))};
  param.vn1={(int)round(profile.Get(31)),(int)round(profile.Get(32)),(int)round(profile.Get(33)),(int)round(profile.Get(38))};

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
  param.nM0=round(profile.Get(9));

  param.beta_sum0=profile.Get(34);
  param.beta_pow0=profile.Get(35);
  param.beta_add0=profile.Get(36);

  param.beta_sum1=profile.Get(34);
  param.beta_pow1=profile.Get(35);
  param.beta_add1=profile.Get(36);

  param.bias_mu0=profile.Get(43);
  param.bias_mu1=profile.Get(44);

  param.bias_scale0=param.bias_scale1=std::round(profile.Get(45));
  //param.mix_wmax=profile.Get(41);
  //param.bias_scale1=std::round(profile.Get(41));

  //param.nM0 = std::min(std::max(0,param.nB-param.nS1),param.nM0);

  //if (param.nS1==0) param.nS1=1;
  param.ch_ref=0;
  if (param.nS1 < 0) {
    param.nS1 = -param.nS1;
    param.ch_ref=1;
  };// else if (param.nS1==0) param.nS1=1;
}

void FrameCoder::PredictFrame(const SacProfile &profile,tch_samples &error,int from,int numsamples,bool optimize)
{

  Predictor::tparam param;
  SetParam(param,profile,optimize);
  Predictor pr(param);

  auto eprocess=[&](int ch_p,int ch,int32_t val,int idx) {
      double pd=pr.predict(ch_p);
      int32_t pi=std::clamp((int32_t)std::round(pd),framestats[ch].minval,framestats[ch].maxval);
      if (!optimize) pred[ch][idx]=pi+framestats[ch].mean;
      error[ch][idx]=val-pi; // needed for cost-function within optimize
      pr.update(ch_p,val);
  };


  if (numchannels_==1) {
    const auto *src=&samples[0][from];
    for (int idx=0;idx<numsamples;idx++)
    {
      pr.fillbuf_ch0(src,idx,src,idx);
      eprocess(0,0,src[idx],idx);
    }
  } else if (numchannels_==2) {
    int ch0=param.ch_ref;
    int ch1=1-ch0;

    const auto *src0=&samples[ch0][from];
    const auto *src1=&samples[ch1][from];

    int idx0=0,idx1=0;
    while (idx0<numsamples || idx1<numsamples)
    {
      if (idx0<numsamples) {
        pr.fillbuf_ch0(src0,idx0,src1,idx1);
        eprocess(0,ch0,src0[idx0],idx0);
        idx0++;
      }
      if (idx0>=param.nS1) {
        pr.fillbuf_ch1(src0,src1,idx1,numsamples);
        eprocess(1,ch1,src1[idx1],idx1);
        idx1++;
      }
    }
  }
}

void FrameCoder::UnpredictFrame(const SacProfile &profile,int numsamples)
{
  Predictor::tparam param;
  SetParam(param,profile,false);
  Predictor pr(param);

  auto dprocess=[&](int ch_p,int ch,int32_t *dst,int idx) {
    const double pd=pr.predict(ch_p);
    const int32_t pi=std::clamp((int32_t)round(pd),framestats[ch].minval,framestats[ch].maxval);


    if (framestats[ch].enc_mapped)
      dst[idx]=pi+framestats[ch].mymap.Unmap(pi+framestats[ch].mean,error[ch][idx]);
    else
      dst[idx]=pi+error[ch][idx];

    pr.update(ch_p,dst[idx]);
  };

  if (numchannels_==1) {
    auto *dst=&samples[0][0];
    for (int idx=0;idx<numsamples;idx++)
    {
      pr.fillbuf_ch0(dst,idx,dst,idx);
      dprocess(0,0,dst,idx);
    }
  } else if (numchannels_==2) {
    int ch0=param.ch_ref;
    int ch1=1-ch0;

    auto *dst0=&samples[ch0][0];
    auto *dst1=&samples[ch1][0];
    int idx0=0,idx1=0;
    while (idx0<numsamples || idx1<numsamples)
    {
      if (idx0<numsamples) {
        pr.fillbuf_ch0(dst0,idx0,dst1,idx1);
        dprocess(0,ch0,dst0,idx0);
        idx0++;
      }
      if (idx0>=param.nS1) {
        pr.fillbuf_ch1(dst0,dst1,idx1,numsamples);
        dprocess(1,ch1,dst1,idx1);
        idx1++;
      }
    }
  }

  // add mean
  for (int ch=0;ch<numchannels_;ch++) {
    if (framestats[ch].mean!=0)
      for (int i=0;i<numsamples;i++) samples[ch][i]+=framestats[ch].mean;
  }
}

int FrameCoder::EncodeMonoFrame_Normal(int ch,int numsamples,BufIO &buf)
{
  buf.Reset();
  RangeCoderSH rc(buf);
  rc.Init();

  BitplaneCoder bc(framestats[ch].maxbpn,numsamples);
  int32_t *psrc=&(s2u_error[ch][0]);
  bc.Encode(rc.encode_p1,psrc);
  rc.Stop();
  return buf.GetBufPos();
}

int FrameCoder::EncodeMonoFrame_Mapped(int ch,int numsamples,BufIO &buf)
{
  buf.Reset();

  RangeCoderSH rc(buf);
  rc.Init();

  BitplaneCoder bc(framestats[ch].maxbpn_map,numsamples);

  MapEncoder me(rc,framestats[ch].mymap.usedl,framestats[ch].mymap.usedh);
  me.Encode();
  bc.Encode(rc.encode_p1,&(s2u_error_map[ch][0]));
  rc.Stop();
  return buf.GetBufPos();
}

double FrameCoder::CalcRemapError(int ch, int numsamples)
{
    std::vector<int32_t>emap(numsamples);
    int32_t emax_map=1;
    for (int i=0;i<numsamples;i++) {
      int32_t map_e=framestats[ch].mymap.Map(pred[ch][i],error[ch][i]);
      int32_t map_ue=MathUtils::S2U(map_e);
      emap[i]=map_e;
      s2u_error_map[ch][i]=map_ue;
      if (map_ue>emax_map) emax_map=map_ue;
    }
    framestats[ch].maxbpn_map=ilogb(emax_map);

    CostL1 cost;

    double ent1 = cost.Calc(std::span{&error[ch][0],static_cast<unsigned>(numsamples)});
    double ent2 = cost.Calc(std::span{&emap[0],static_cast<unsigned>(numsamples)});
    double r=1.0;
    if (ent2!=0.0) r=ent1/ent2;
    if (opt.verbose_level>0) std::cout << "  cost pcm-model: " << ent1 << ' ' << ent2 << ' ' << r << '\n';
    return r;
}

void FrameCoder::EncodeMonoFrame(int ch,int numsamples)
{
  if (opt.sparse_pcm==0) {
    EncodeMonoFrame_Normal(ch,numsamples,enc_temp1[ch]);
    framestats[ch].enc_mapped=false;
    encoded[ch]=enc_temp1[ch];
  } else {
    double r = CalcRemapError(ch,numsamples);
    int size_normal=EncodeMonoFrame_Normal(ch,numsamples,enc_temp1[ch]);
    framestats[ch].enc_mapped=false;
    encoded[ch]=enc_temp1[ch];

    if (r > 1.05)
    {
      int size_mapped=EncodeMonoFrame_Mapped(ch,numsamples,enc_temp2[ch]);
      if (size_mapped<size_normal)
      {
        if (opt.verbose_level>0) {
          std::cout << "  sparse frame " << size_normal << " -> " << size_mapped << " (" << (size_mapped-size_normal) << ")\n";
        }
        framestats[ch].enc_mapped=true;
        encoded[ch]=enc_temp2[ch];
      }
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

  BitplaneCoder bc(framestats[ch].maxbpn,numsamples);
  bc.Decode(rc.decode_p1,dst);
  rc.Stop();
}


void FrameCoder::PrintProfile(SacProfile &profile)
{
    Predictor::tparam param;
    SetParam(param,profile);

    std::cout << '\n';
    std::cout << "lpc (nA " << std::round(profile.Get(24)) << " nM0 " << std::round(profile.Get(9));
    std::cout << ") (nB " << std::round(profile.Get(25)) << " nS0 " << std::round(profile.Get(26)) << " nS1 " << std::round(profile.Get(27)) << ")\n";
    std::cout << "lpc nu " << param.ols_nu0 << ' ' << param.ols_nu1 << '\n';
    std::cout << "lpc cov0 " << param.beta_sum0 << ' ' << param.beta_pow0 << ' ' << param.beta_add0 << "\n";
    std::cout << "lms0 ";
    for (int i=28;i<=30;i++) std::cout << std::round(profile.Get(i)) << ' ';
    std::cout << std::round(profile.Get(37));
    std::cout << '\n';
    std::cout << "lms1 ";
    for (int i=31;i<=33;i++) std::cout << std::round(profile.Get(i)) << ' ';
    std::cout << std::round(profile.Get(38));
    std::cout << '\n';
    std::cout << "mu ";
    for (std::size_t i=0;i<std::size(param.vmu0);i++)
      std::cout << (param.vmu0[i]*param.vn0[i]) << ' ';
    std::cout << '\n';
    std::cout << "mu_decay ";
    for (auto &x : param.vmudecay0)
      std::cout << x << ' ';
    std::cout << '\n';
    std::cout << "pow_decay ";
    for (const auto &x : param.vpowdecay0)
      std::cout << x << ' ';
    std::cout << '\n';

    std::cout << "mu mix beta " << param.mu_mix_beta0 << " " << param.mu_mix_beta1 << '\n';
    std::cout << "ch-ref " << param.ch_ref << "\n";
    std::cout << "bias mu " << param.bias_mu0 << ", " << param.bias_mu1 << " scale " << (1<<param.bias_scale0) << ' ' << (1<<param.bias_scale1) << '\n';
}

double FrameCoder::GetCost(const CostFunction *func,const tch_samples &samples,std::size_t samples_to_optimize) const
{
  // return a span over samples
  const auto span_ch = [=,this](int ch){
    return std::span{&samples[ch][0],samples_to_optimize};
  };

  double cost=0.0;
  if (opt.mt_mode>1 && numchannels_>1) {

    std::vector <std::future<double>> threads;
    for (int ch=0;ch<numchannels_;ch++)
        threads.emplace_back(std::async([=]{return func->Calc(span_ch(ch));}));

    for (auto &thread : threads)
      cost += thread.get();

  } else {
    for (int ch=0;ch<numchannels_;ch++)
      cost += func->Calc(span_ch(ch));
  }
  return cost;
}

void FrameCoder::Optimize(const FrameCoder::toptim_cfg &ocfg,SacProfile &profile,const std::vector<int>&params_to_optimize)
{
  int samples_to_optimize=std::min(numsamples_,static_cast<int>(std::ceil(framesize_*ocfg.fraction)));
  const int start_pos=(numsamples_-samples_to_optimize)/2;

  CostFunction *CostFunc=nullptr;
  switch (ocfg.optimize_cost)  {
    case FrameCoder::SearchCost::L1:CostFunc=new CostL1();break;
    case FrameCoder::SearchCost::RMS:CostFunc=new CostRMS();break;
    case FrameCoder::SearchCost::Golomb:CostFunc=new CostGolomb();break;
    case FrameCoder::SearchCost::Entropy:CostFunc=new CostEntropy();break;
    case FrameCoder::SearchCost::Bitplane:CostFunc=new CostBitplane();break;
    default:std::cerr << "  error: unknown FramerCoder::CostFunction\n";return;
  }

  const int ndim=params_to_optimize.size();
  vec1D xstart(ndim); // starting vector
  Opt::box_const pb(ndim); // set constraints
  for (int i=0;i<ndim;i++) {
    pb[i].xmin=profile.coefs[params_to_optimize[i]].vmin;
    pb[i].xmax=profile.coefs[params_to_optimize[i]].vmax;
    xstart[i]=profile.coefs[params_to_optimize[i]].vdef;
  }

  auto cost_func=[&](const vec1D &x) {
    // create thread safe copies for error and profile
    tch_samples tmp_error(numchannels_,std::vector<int32_t>(samples_to_optimize));
    SacProfile tmp_profile=profile;

    for (int i=0;i<ndim;i++) tmp_profile.coefs[params_to_optimize[i]].vdef=x[i];

    PredictFrame(tmp_profile,tmp_error,start_pos,samples_to_optimize,true);
    return GetCost(CostFunc,tmp_error,samples_to_optimize);
  };

  if (opt.verbose_level>0) {
    std::string opt_str="DDS";
    if (ocfg.optimize_search==FrameCoder::SearchMethod::DE) opt_str="DE";
    std::cout << "\n " << opt_str << " " << ocfg.maxnfunc << "= ";
  }

  std::unique_ptr<Opt> myOpt;

  if (ocfg.optimize_search==FrameCoder::SearchMethod::DDS)
    myOpt = std::make_unique<OptDDS>(ocfg.dds_cfg,pb,opt.verbose_level);
  else if (ocfg.optimize_search==FrameCoder::SearchMethod::DE)
    myOpt = std::make_unique<OptDE>(ocfg.de_cfg,pb,opt.verbose_level);

  Opt::ppoint ret = myOpt->run(cost_func,xstart);

  // save optimal vector to baseprofile
  for (int i=0;i<ndim;i++)
    profile.coefs[params_to_optimize[i]].vdef=ret.second[i];

  if (opt.verbose_level>0) {
    PrintProfile(profile);
  }

  delete CostFunc;
}

void FrameCoder::CnvError_S2U(tch_samples &error,int numsamples)
{
  for (int ch=0;ch<numchannels_;ch++)
  {
    int32_t emax=1;
    for (int i=0;i<numsamples;i++) {
      const int32_t e_s2u=MathUtils::S2U(error[ch][i]);
      if (e_s2u>emax) emax=e_s2u;
      s2u_error[ch][i]=e_s2u;
    }
    framestats[ch].maxbpn=ilogb(emax);
  }
}

void FrameCoder::Predict()
{
  for (int ch=0;ch<numchannels_;ch++)
  {
    AnalyseMonoChannel(ch,numsamples_);
    if (opt.sparse_pcm) {
      framestats[ch].mymap.Reset();
      framestats[ch].mymap.Analyse(&(samples[ch][0]),numsamples_);
    }
    if (opt.zero_mean==0) {
      framestats[ch].mean = 0;
    } else if (framestats[ch].mean!=0) {
      for (int i=0;i<numsamples_;i++) samples[ch][i] -= framestats[ch].mean;
      framestats[ch].minval -= framestats[ch].mean;
      framestats[ch].maxval -= framestats[ch].mean;
    }
  }

  if (opt.optimize)
  {
    // reset profile params
    // otherwise: starting point for optimization is the best point from the last frame
    if (opt.ocfg.reset)
      base_profile.LoadBaseProfile();

    // optimize all params
    std::vector<int>lparam_base(base_profile.coefs.size());
    std::iota(std::begin(lparam_base),std::end(lparam_base),0);

    Optimize(opt.ocfg,base_profile,lparam_base);
  }
  PredictFrame(base_profile,error,0,numsamples_,false);
  CnvError_S2U(error,numsamples_);
}

void FrameCoder::Unpredict()
{
  UnpredictFrame(base_profile,numsamples_);
}

void FrameCoder::Encode()
{
  if (opt.mt_mode && numchannels_>1)  {
    std::vector <std::jthread> threads;
    for (int ch=0;ch<numchannels_;ch++)
      threads.emplace_back(std::jthread(&FrameCoder::EncodeMonoFrame,this,ch,numsamples_));
  } else
    for (int ch=0;ch<numchannels_;ch++) EncodeMonoFrame(ch,numsamples_);
}

void FrameCoder::Decode()
{
  if (opt.mt_mode && numchannels_>1) {
    std::vector <std::jthread> threads;
    for (int ch=0;ch<numchannels_;ch++)
      threads.emplace_back(std::jthread(&FrameCoder::DecodeMonoFrame,this,ch,numsamples_));
  } else
    for (int ch=0;ch<numchannels_;ch++)
      DecodeMonoFrame(ch,numsamples_);
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

int FrameCoder::WriteBlockHeader(std::fstream &file, const std::vector<SacProfile::FrameStats> &framestats,int ch)
{
  uint8_t buf[32];
  BitUtils::put32LH(buf,framestats[ch].blocksize);
  BitUtils::put32LH(buf+4,static_cast<uint32_t>(framestats[ch].mean));
  BitUtils::put32LH(buf+8,static_cast<uint32_t>(framestats[ch].minval));
  BitUtils::put32LH(buf+12,static_cast<uint32_t>(framestats[ch].maxval));
  uint16_t flag=0;
  if (framestats[ch].enc_mapped) {
    flag|=(1<<9);
    flag|=framestats[ch].maxbpn_map;
  } else {
    flag|=framestats[ch].maxbpn;
  }
  BitUtils::put16LH(buf+16,flag);
  file.write(reinterpret_cast<char*>(buf),18);
  return 18;
}

int FrameCoder::ReadBlockHeader(std::fstream &file, std::vector<SacProfile::FrameStats> &framestats,int ch)
{
  uint8_t buf[32];
  file.read(reinterpret_cast<char*>(buf),18);

  framestats[ch].blocksize=BitUtils::get32LH(buf);
  framestats[ch].mean=static_cast<int32_t>(BitUtils::get32LH(buf+4));
  framestats[ch].minval=static_cast<int32_t>(BitUtils::get32LH(buf+8));
  framestats[ch].maxval=static_cast<int32_t>(BitUtils::get32LH(buf+12));
  uint16_t flag=BitUtils::get16LH(buf+16);
  if (flag>>9) framestats[ch].enc_mapped=true;
  else framestats[ch].enc_mapped=false;
  framestats[ch].maxbpn=flag&0xff;
  return 18;
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

  if (numsamples) {
    int64_t sum=0;
    for (int i=0;i<numsamples;i++) {
        sum += src[i];
    }
    framestats[ch].mean = (int)std::floor(sum / (double)numsamples);

    int32_t minval = std::numeric_limits<int32_t>::max();
    int32_t maxval = std::numeric_limits<int32_t>::min();
    for (int i=0;i<numsamples;i++) {
      const int32_t val=src[i];
      if (val>maxval) maxval=val;
      if (val<minval) minval=val;
    }
    framestats[ch].minval = minval;
    framestats[ch].maxval = maxval;
    if (opt.verbose_level>0) {
      std::cout << "  ch" << ch << " samples=" << numsamples;
      std::cout << ",mean=" << framestats[ch].mean << ",min=" << framestats[ch].minval << ",max=" << framestats[ch].maxval << "\n";
    }
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
  profile_tmp.LoadBaseProfile();
  const int size_profile_bytes=profile_tmp.coefs.size()*4;

  int frame_num=1;
  int coef_hdr_size=0;
  int block_hdr_size=0;
  while (mySac.file.tellg()<fsize) {
    uint8_t buf[12];
    mySac.file.read(reinterpret_cast<char*>(buf),4);
    int numsamples=BitUtils::get32LH(buf);
    std::cout << "Frame " << frame_num << ": " << numsamples << " samples "<< std::endl;

    mySac.file.seekg(size_profile_bytes,std::ios_base::cur); // skip profile coefs
    coef_hdr_size += size_profile_bytes;


    for (int ch=0;ch<mySac.getNumChannels();ch++) {
      int num_bytes=FrameCoder::ReadBlockHeader(mySac.file, framestats, ch);
      block_hdr_size += num_bytes;
      std::cout << "  Channel " << ch << ": " << framestats[ch].blocksize << " bytes\n";
      std::cout << "    Bpn: " << framestats[ch].maxbpn << ", sparse_pcm: " << (framestats[ch].enc_mapped) << std::endl;
      std::cout << "    mean: " << framestats[ch].mean << ", min: " << framestats[ch].minval << ", max: " << framestats[ch].maxval << std::endl;
      mySac.file.seekg(framestats[ch].blocksize, std::ios_base::cur);
    }
    frame_num++;
  }
  std::cout << "Frames   " << (frame_num-1) << '\n';
  std::cout << "Hdr_size " << (coef_hdr_size+block_hdr_size) << " (coefs " << coef_hdr_size << ",block " << block_hdr_size << ")\n";
}


std::pair<double,double> Codec::AnalyseSparse(std::span<const int32_t> buf)
{
  SparsePCM spcm;
  spcm.Analyse(buf);

  return {spcm.fraction_used,spcm.fraction_cost};
}

void Codec::PushState(std::vector<Codec::tsub_frame> &sub_frames,Codec::tsub_frame &curframe,int min_frame_length,int block_state=-1,int samples_block=0)
{
  if (block_state==curframe.state)
    curframe.length+=samples_block;
  else {
    if (curframe.length < min_frame_length && sub_frames.size()) // extend
    {
      sub_frames.back().length+=curframe.length;
    } else {
      if (opt_.verbose_level>1)
        std::cout << "push subframe of length " << curframe.length << " samples\n";
      sub_frames.push_back(curframe);

      if (samples_block) {
        curframe.state=block_state; // set new blockstate
        curframe.start+=curframe.length;
        curframe.length=samples_block;
      }
    }
  }
}

std::vector<Codec::tsub_frame> Codec::Analyse(const std::vector <std::vector<int32_t>>&samples,int blocksamples,int min_frame_length,int samples_read)
{
  std::vector<Codec::tsub_frame> sub_frames;

  int samples_processed=0;
  int nblock=0;

  Codec::tsub_frame curframe;

  while (samples_processed < samples_read)
  {
    int samples_left = samples_read-samples_processed;
    int samples_block = std::min(blocksamples,samples_left);
    double avg_cost=0,avg_used=0;
    for (unsigned ch=0;ch<samples.size();ch++)
    {
      auto [fused,fcost]=AnalyseSparse(std::span{&samples[ch][samples_processed],static_cast<unsigned>(samples_block)});
      avg_cost+=fcost;
      avg_used+=fused;
    }
    avg_cost /= (double)samples.size();
    avg_used /= (double)samples.size();
    int block_state=(avg_cost>1.35); // high threshold
    if (opt_.verbose_level>1) {
      std::cout << "  analyse block " << nblock << ' ' << samples_block << " sparse " << block_state << " (" << avg_cost << "," << avg_used << ")\n";
    }

    if (nblock==0)
    {
      curframe.state = block_state;
      curframe.length=samples_block;
      curframe.start=0;
    } else
      PushState(sub_frames,curframe,min_frame_length,block_state,samples_block);

    samples_processed += samples_block;
    nblock++;
  }

  if (curframe.length)
    PushState(sub_frames,curframe,min_frame_length);

  if (samples_processed != samples_read)
    std::cerr << "  warning: samples_processed != samples_read (" << samples_processed << "," << samples_read << ")\n";

  if (opt_.verbose_level>1) std::cout << "sub_frames\n";
  int64_t nlen=0;
  for (const auto &frame : sub_frames) {
    if (opt_.verbose_level>1) std::cout << "  " << frame.start << ' ' << frame.length << ' ' << frame.state << '\n';
    nlen+=frame.length;
  }
  if (nlen!=samples_read)
    std::cerr << "  warning: nlen != samples_read\n";
  return sub_frames;
}

void Codec::EncodeFile(Wav &myWav,Sac &mySac)
{
  uint32_t max_framesize=static_cast<uint32_t>(opt_.max_framelen)*myWav.getSampleRate();

  const int numchannels=myWav.getNumChannels();

  FrameCoder myFrame(numchannels,max_framesize,opt_);

  mySac.mcfg.max_framelen = opt_.max_framelen;

  mySac.WriteSACHeader(myWav);
  std::streampos hdrpos = mySac.file.tellg();
  mySac.WriteMD5(myWav.md5ctx.digest);
  myWav.InitFileBuf(max_framesize);

  Timer gtimer,ltimer;
  double time_prd=0,time_enc=0;

  gtimer.start();
  int samplescoded=0;
  int samplestocode=myWav.getNumSamples();
  std::vector<std::vector<int32_t>> csamples(myWav.getNumChannels(),std::vector<int32_t>(max_framesize));

  while (samplestocode>0) {
      int samplesread=myWav.ReadSamples(csamples,max_framesize);

      std::vector<Codec::tsub_frame> sub_frames;
      if (opt_.adapt_block) {
        int block_len=myWav.getSampleRate()*3;
        int min_frame_len=myWav.getSampleRate()*3;
        sub_frames=Analyse(csamples,block_len,min_frame_len,samplesread);
      } else {
        sub_frames.push_back({0,0,samplesread});
      }

      for (auto &subframe:sub_frames)
      {
        if (opt_.verbose_level)
          std::cout << "frame " << subframe.start << " state " << subframe.state << " len " << subframe.length << '\n';

        for (int ch=0;ch<myWav.getNumChannels();ch++)
          std::copy_n(&csamples[ch][subframe.start],subframe.length,&myFrame.samples[ch][0]);

        myFrame.SetNumSamples(subframe.length);

        ltimer.start();myFrame.Predict();ltimer.stop();time_prd+=ltimer.elapsedS();
        ltimer.start();myFrame.Encode();ltimer.stop();time_enc+=ltimer.elapsedS();
        myFrame.WriteEncoded(mySac);

        samplescoded+=subframe.length;
        PrintProgress(samplescoded,myWav.getNumSamples());
        samplestocode-=subframe.length;
      }
  }
  MD5::Finalize(&myWav.md5ctx);
  gtimer.stop();
  double time_total=gtimer.elapsedS();
  if (time_total>0.)   {
     double rprd=time_prd*100./time_total;
     double renc=time_enc*100./time_total;
     std::cout << "\n  Timing:  pred " << miscUtils::ConvertFixed(rprd,2) << "%, ";
     std::cout << "enc " << miscUtils::ConvertFixed(renc,2) << "%, ";
     std::cout << "misc " << miscUtils::ConvertFixed(100.-rprd-renc,2) << "%" << std::endl;
  }
  std::cout << "  MD5:     ";
  for (auto x : myWav.md5ctx.digest) std::cout << std::hex << (int)x;
  std::cout << std::dec << '\n';

  std::streampos eofpos = mySac.file.tellg();
  mySac.file.seekg(hdrpos);
  mySac.WriteMD5(myWav.md5ctx.digest);
  mySac.file.seekg(eofpos);
}

void Codec::DecodeFile(Sac &mySac,Wav &myWav)
{
  const Sac::sac_cfg &cfg=mySac.mcfg;
  myWav.InitFileBuf(cfg.max_framesize);
  mySac.UnpackMetaData(myWav);
  myWav.WriteHeader();

  opt_.max_framelen=cfg.max_framelen;
  FrameCoder myFrame(mySac.getNumChannels(),cfg.max_framesize,opt_);

  int64_t data_nbytes=0;
  int samplestodecode=mySac.getNumSamples();
  int samplesdecoded=0;
  while (samplestodecode>0) {
    myFrame.ReadEncoded(mySac);
    myFrame.Decode();
    myFrame.Unpredict();
    data_nbytes += myWav.WriteSamples(myFrame.samples,myFrame.GetNumSamples());

    samplesdecoded+=myFrame.GetNumSamples();
    PrintProgress(samplesdecoded,myWav.getNumSamples());
    samplestodecode-=myFrame.GetNumSamples();
  }
  // pad odd sized data chunk
  if (data_nbytes&1) myWav.WriteData(std::vector<uint8_t>{0},1);
  myWav.WriteHeader();
}
