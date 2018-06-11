#include "libsac.h"
#include "../common/timer.h"

FrameCoder::FrameCoder(int numchannels,int framesize,int profile,int sparse_pcm)
:numchannels_(numchannels),framesize_(framesize),sparse_pcm_(sparse_pcm)
{
  if (profile==0)   {
    baseprofile.Init(2,0);
    baseprofile.Set(0,0.99,0.999,0.996);
    baseprofile.Set(1,0.1,2.0,0.5);
  } else if (profile==1) {
    baseprofile.Init(6,1);
    baseprofile.Set(0,0.99,0.9999,0.998);
    baseprofile.Set(1,0.000001,0.001,0.0004);
    baseprofile.Set(2,0.000001,0.001,0.0001);
    baseprofile.Set(3,0.00001,0.01,0.003);
    baseprofile.Set(4,0.0001,0.1 ,0.01);
    baseprofile.Set(5,0.99,0.9999,0.998);
  } else if (profile==2) {
    baseprofile.Init(8,2);
    baseprofile.Set(0,0.98,0.9999,0.997);
    baseprofile.Set(1,0.0000001,0.0001,0.00001);
    baseprofile.Set(2,0.0000001,0.0001,0.000005);

    baseprofile.Set(3,0.000001,0.001,0.0005);
    baseprofile.Set(4,0.000001,0.001,0.0001);
    baseprofile.Set(5,0.0001,0.01 ,0.003);
    baseprofile.Set(6,0.001,0.1   ,0.01);
    baseprofile.Set(7,0.99,0.9999,0.998);
  }
  profile_size_bytes_=baseprofile.coefs.size()*4;

  framestats.resize(numchannels);
  samples.resize(numchannels);
  err0.resize(numchannels);
  err1.resize(numchannels);
  error.resize(numchannels);
  s2u_error.resize(numchannels);
  s2u_error_map.resize(numchannels);
  for (int i=0;i<numchannels;i++) {
    samples[i].resize(framesize);
    err0[i].resize(framesize);
    err1[i].resize(framesize);
    error[i].resize(framesize);
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
    /*int tshift=0;
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
    }*/
    //cout << "total shift: " << tshift << endl;

    mean=mean/numsamples;
    /*if (mean) {
       minval-=mean;
       maxval-=mean;
       for (int i=0;i<numsamples;i++) src[i]-=mean;
    }*/
    framestats[ch].minval=minval;
    framestats[ch].maxval=maxval;
    framestats[ch].mean=mean;
    std::cout << "mean: " << mean << ", min: " << minval << " ,max: " << maxval << std::endl;
  }
}

StereoPredictor *FrameCoder::AllocProfile(const SacProfile &profile)
{
  StereoPredictor *pred=nullptr;
  switch (profile.type) {
    case 0:pred=new StereoFast(profile);break;
    case 1:pred=new StereoNormal(profile);break;
    case 2:pred=new StereoHigh(profile);break;
    default:break;
  }
  return pred;
}

int Signed2Unsigned(int32_t pred,int32_t error,int32_t minval,int32_t maxval)
{
  int32_t abserr=std::abs(error);
  if (abserr) {
    if (pred <= 0) {
      int32_t diff=(-minval)+pred;
      if (error<0) return (2*(-error));
      else {
        if (error<=diff) return (2*error-1);
        else return (diff+error);
      }
    } else {
      int32_t diff=maxval-pred;
      if (error>0) return(2*error-1);
      else {
        if (abserr<=diff) return (2*abserr);
        else return (diff+abserr);
      }
    }
  } return 0;
}

void FrameCoder::PredictStereoFrame(const SacProfile &profile,int ch0,int ch1,int from,int numsamples)
{
 StereoPredictor *pred=AllocProfile(profile);

 if (pred!=nullptr) {
  const int32_t *src0=&(samples[ch0][from]);
  const int32_t *src1=&(samples[ch1][from]);
  int32_t *dst0=&(error[ch0][0]);
  int32_t *dst1=&(error[ch1][0]);
  int32_t emax0=0,emax1=0;
  int32_t emax0_map=0,emax1_map=0;

  for (int i=0;i<numsamples;i++) {
    int32_t p0=pred->Predict();
    pred->Update(src0[i]);
    int32_t p1=pred->Predict();
    pred->Update(src1[i]);
    p0=clamp(p0,-(1<<15),1<<15);
    p1=clamp(p1,-(1<<15),1<<15);

    dst0[i]=src0[i]-p0;
    dst1[i]=src1[i]-p1;

    int map_a=framestats[0].mymap.Map(p0,dst0[i]); // remap the error
    int map_b=framestats[1].mymap.Map(p1,dst1[i]);

    #if 0 //debug code
      int unmap_a=framestats[0].mymap.Unmap(pred0,map_a);
      int unmap_b=framestats[1].mymap.Unmap(pred1,map_b);
      if (unmap_a != dst0[i]) cout << unmap_a << " " << dst0[i] << endl;
      if (unmap_b != dst1[i]) cout << unmap_b << " " << dst1[i] << endl;
    #endif

    int32_t e0=MathUtils::S2U(dst0[i]);
    int32_t e1=MathUtils::S2U(dst1[i]);

    //int e0x=Signed2Unsigned(p0,dst0[i],-(1<<15),1<<15);
    //int e1x=Signed2Unsigned(p1,dst1[i],-(1<<15),1<<15);

    //if (e0x < e0) std::cout << e0 << " " << e0x << std::endl;

    int32_t e0_map=MathUtils::S2U(map_a);
    int32_t e1_map=MathUtils::S2U(map_b);

    s2u_error[ch0][i]=e0;
    s2u_error[ch1][i]=e1;

    s2u_error_map[ch0][i]=e0_map;
    s2u_error_map[ch1][i]=e1_map;

    if (e0>emax0) emax0=e0;
    if (e1>emax1) emax1=e1;

    if (e0_map>emax0_map) emax0_map=e0_map;
    if (e1_map>emax1_map) emax1_map=e1_map;
  }

  framestats[0].maxbpn=MathUtils::iLog2(emax0);
  framestats[1].maxbpn=MathUtils::iLog2(emax1);

  framestats[0].maxbpn_map=MathUtils::iLog2(emax0_map);
  framestats[1].maxbpn_map=MathUtils::iLog2(emax1_map);
  delete pred;
 }
}

void FrameCoder::UnpredictStereoFrame(const SacProfile &profile,int ch0,int ch1,int numsamples)
{
  StereoPredictor *pred=AllocProfile(profile);

  const int32_t *src0=&(error[ch0][0]);
  const int32_t *src1=&(error[ch1][0]);
  int32_t *dst0=&(samples[ch0][0]);
  int32_t *dst1=&(samples[ch1][0]);

  for (int i=0;i<numsamples;i++) {
    int32_t p0=pred->Predict();
    p0=clamp(p0,-(1<<15),1<<15);

    if (framestats[0].enc_mapped) dst0[i]=p0+framestats[0].mymap.Unmap(p0,src0[i]);
    else dst0[i]=p0+src0[i];
    pred->Update(dst0[i]);

    int32_t p1=pred->Predict();
    p1=clamp(p1,-(1<<15),1<<15);

    if (framestats[1].enc_mapped) dst1[i]=p1+framestats[1].mymap.Unmap(p1,src1[i]);
    else dst1[i]=p1+src1[i];
    pred->Update(dst1[i]);
  }

  delete pred;
}

int FrameCoder::EncodeMonoFrame_Normal(int ch,int numsamples,BufIO &buf)
{
  buf.Reset();
  RangeCoderSH rc(buf);
  rc.Init();
  BitplaneCoder bc(rc,framestats[ch].maxbpn,numsamples);
  int32_t *srca=&(s2u_error[ch][0]);
  bc.Encode(srca);
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
  //cout << "mapsize: " << buf.GetBufPos() << " Bytes\n";

  bc.Encode(&(s2u_error_map[ch][0]));
  rc.Stop();
  return buf.GetBufPos();
}


void FrameCoder::EncodeMonoFrame(int ch,int numsamples)
{
  if (sparse_pcm_==0) {
    int size_normal=EncodeMonoFrame_Normal(ch,numsamples,enc_temp1[ch]);
    framestats[ch].enc_mapped=false;
    encoded[ch]=enc_temp1[ch];
  } else {
    int size_normal=EncodeMonoFrame_Normal(ch,numsamples,enc_temp1[ch]);
    int size_mapped=EncodeMonoFrame_Mapped(ch,numsamples,enc_temp2[ch]);
    if (size_normal<size_mapped)
    {
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
    std::cout << buf.GetBufPos() << std::endl;
  }

  BitplaneCoder bc(rc,framestats[ch].maxbpn,numsamples);
  bc.Decode(dst);
  rc.Stop();
}

class CostMeanRMS : public CostFunction {
  public:
      double Calc(const int32_t *buf,int numsamples)
      {
        if (numsamples) {
          double pow=0.0;
          for (int i=0;i<numsamples;i++) pow+=std::abs(buf[i]);
          return pow/static_cast<double>(numsamples);
        } else return 0.;
      }
};

// estimate number of needed bits with a simple golomb model
class CostGolomb : public CostFunction {
  public:
      CostGolomb():mean_err(0.98){};
      double Calc(const int32_t *buf,int numsamples)
      {
        double nbits=0;
        if (numsamples) {
          for (int i=0;i<numsamples;i++) {
            int32_t val=MathUtils::S2U(buf[i]);
            int m=std::max(static_cast<int>(mean_err.Get()),1);
            int q=val/m;
            //int r=val-q*m;
            nbits+=(q+1);
            if (m>1) {
              int b=ceil(log(m)/log(2));
              nbits+=b;
            }
            mean_err.Update(val);
          }
          return nbits/(double)numsamples;
        } else return 0;
      }
  private:
    RunExp mean_err;
};

double FrameCoder::GetCost(SacProfile &profile,CostFunction *func,int coef,double testval,int start_sample,int samples_to_optimize)
{
  const int32_t *err0=&(error[0][0]);
  const int32_t *err1=&(error[1][0]);

  SacProfile testprofile=profile;
  testprofile.coefs[coef].vdef=testval;
  PredictStereoFrame(testprofile,0,1,start_sample,samples_to_optimize);
  double c1=func->Calc(err0,samples_to_optimize);
  double c2=func->Calc(err1,samples_to_optimize);
  return (c1+c2)/2.;
}

int SolveQuadratic(double x0,double x1,double x2,double y0,double y1,double y2,double &xm)
{
  double det=(x0-x1)*(x1*x1-x2*x2)+(x2-x1)*(x0*x0-x1*x1);
  if (fabs(det)>0.0) {
    double b=((y0-y1)*(x1*x1-x2*x2)+(y2-y1)*(x0*x0-x1*x1))/det;
    double a=(y0-y1-b*(x0-x1))/(x0*x0-x1*x1);
    double c=y0-a*(x0*x0)-b*x0;

    if (a>0.) {xm=-b/(2*a);return 0;};
  }
  return 1;
}

void FrameCoder::Optimize(SacProfile &profile)
{
  int blocksize=4;
  int samples_to_optimize=numsamples_/blocksize;
  //cout << numsamples << "->" << samples_to_optimize << " " << endl;
  int start=(numsamples_-samples_to_optimize)/2;
  //cout << start << "->" << start+samples_to_optimize << " " << endl;
  const int32_t *err0=&(error[0][0]);
  const int32_t *err1=&(error[1][0]);

  CostGolomb CostFunc;

  #if 0
  double tx[4],ty[4];
  for (int coef=0;coef<2;coef++) {
    tx[0]=profile.coefs[coef].vdef;
    tx[1]=profile.coefs[coef].vmin;
    tx[2]=profile.coefs[coef].vmax;
    tx[3]=0.;
    ty[0]=GetCost(profile,&CostFunc,coef,tx[0],start,samples_to_optimize);
    ty[1]=GetCost(profile,&CostFunc,coef,tx[1],start,samples_to_optimize);
    ty[2]=GetCost(profile,&CostFunc,coef,tx[2],start,samples_to_optimize);
    ty[3]=std::numeric_limits<double>::max();

    if (SolveQuadratic(tx[0],tx[1],tx[2],ty[0],ty[1],ty[2],tx[3])==0) {
      ty[3]=GetCost(profile,&CostFunc,coef,tx[3],start,samples_to_optimize);
    }

    double xbest=tx[0];
    double ybest=ty[0];
    for (int i=1;i<4;i++)
      if (ty[i]<ybest) {xbest=tx[i];ybest=ty[i];};

    profile.coefs[coef].vdef=xbest;
  }
  #else
  //std::cout << xbest << ":" << ybest << std::endl;

  /*vector <SacProfile::coef>&coefs=testprofile.coefs;
  int n=static_cast<int>(coefs.size());
  vector <double>dx(n),gx(n);
  for (int i=0;i<n;i++) {
      dx[i]=(coefs[i].vmax-coefs[i].vmin)*0.01;
  }*/

  /*double cost,dwcost,tcost,grad;
  testprofile=profile;
  cost=GetCost(testprofile,&CostFunc,0,coefs[0].vdef,start,samples_to_optimize);
  dwcost=GetCost(testprofile,&CostFunc,0,coefs[0].vdef+dx[0],start,samples_to_optimize);
  grad=(dwcost-cost)/dx[0];

  cout << cost << " " << dwcost << endl;
  double mu=0.01;*/

  //cout << (coefs[0].vdef-mu*grad/fabs(grad)) << endl;
  //profile.coefs[0].vdef=x;

  double best_cost=GetCost(profile,&CostFunc,0,profile.coefs[0].vdef,start,samples_to_optimize);

  //for (int k=0;k<2;k++)
  for (int i=0;i<2;i++)
  {
      double best_x=profile.coefs[i].vdef;
      double x=profile.coefs[i].vmin;
      double xinc=(profile.coefs[i].vmax-profile.coefs[i].vmin)/10.;
      //cout << i << ":" << best_x << " -> ";
      while (x<=profile.coefs[i].vmax)
      {
          double cost=GetCost(profile,&CostFunc,i,x,start,samples_to_optimize);
          //cout << x << " " << cost << " " << best_cost << endl;
          if (cost<best_cost) {
            best_x=x;
            best_cost=cost;
          }
          x+=xinc;
      }
      std::cout << best_x << std::endl;
      profile.coefs[i].vdef=best_x;
  }
  #endif
}

void FrameCoder::InterChannel(int ch0,int ch1,int numsamples)
{
  int32_t *src0=&(samples[ch0][0]);
  int32_t *src1=&(samples[ch1][0]);
  uint64_t sum[4],score[4];
  sum[0]=sum[1]=sum[2]=sum[3]=0;
  for (int i=0;i<numsamples;i++) {
    int32_t v0=src0[i];
    int32_t v1=src1[i];

    sum[0]+=abs(v0);
    sum[1]+=abs(v1);
    sum[2]+=abs((v0+v1)>>1);
    sum[3]+=abs(v0-v1);
  }
  score[0]=sum[0]+sum[1];
  score[1]=sum[0]+sum[3];
  score[2]=sum[1]+sum[3];
  score[3]=sum[2]+sum[3];

  if (score[3]*1.5 < score[0]) {
    std::cout << "ms\n";
    for (int i=0;i<numsamples;i++) {
     int32_t v0=src0[i];
     int32_t v1=src1[i];
     src0[i]=(v0+v1)>>1;
     src1[i]=v0-v1;

   };
  }
}

void FrameCoder::Predict()
{
  //for (int ch=0;ch<numchannels_;ch++) AnalyseChannel(ch,numsamples_);
  //AnalyseSpectrum(0,numsamples);
  //for (int ch=0;ch<numchannels;ch++) AnalyseChannel(ch,numsamples);
  //AnalyseSpectrum(1,numsamples);
  for (int ch=0;ch<numchannels_;ch++)
  {
    framestats[ch].mymap.Reset();
    framestats[ch].mymap.Analyse(&(samples[ch][0]),numsamples_);
  }
  //cout << framestats[0].mymap.Compare(framestats[1].mymap) << endl;

  //SacProfile optprofile=baseprofile;
  //Optimize(baseprofile);
  if (numchannels_==2) {
    //InterChannel(0,1,numsamples);
    PredictStereoFrame(baseprofile,0,1,0,numsamples_);
  }// else for (int ch=0;ch<numchannels;ch++) PredictMonoFrame(ch,numsamples);
}

void FrameCoder::Unpredict()
{
  if (numchannels_==2) UnpredictStereoFrame(baseprofile,0,1,numsamples_);
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

void FrameCoder::WriteEncoded(AudioFile &fout)
{
  uint8_t buf[8];
  BitUtils::put32LH(buf,numsamples_);
  fout.file.write(reinterpret_cast<char*>(buf),4);
  std::vector <uint8_t>profile_buf(profile_size_bytes_);
  EncodeProfile(baseprofile,profile_buf);
  fout.file.write(reinterpret_cast<char*>(&profile_buf[0]),profile_size_bytes_);
  //cout << "numsamples: " << numsamples << endl;
  //cout << " filepos:   " << fout.file.tellg() << endl;
  for (int ch=0;ch<numchannels_;ch++) {
    uint32_t blocksize=encoded[ch].GetBufPos();
    BitUtils::put32LH(buf,blocksize);
    uint16_t flag=0;
    if (framestats[ch].enc_mapped) {
       flag|=(1<<9);
       flag|=framestats[ch].maxbpn_map;
    } else {
       flag|=framestats[ch].maxbpn;
    }
    BitUtils::put16LH(buf+4,flag);

    //cout << " blocksize: " << blocksize << endl;
    fout.file.write(reinterpret_cast<char*>(buf),6);
    fout.WriteData(encoded[ch].GetBuf(),blocksize);
  }
}

void FrameCoder::ReadEncoded(AudioFile &fin)
{
  uint8_t buf[8];
  fin.file.read(reinterpret_cast<char*>(buf),4);
  numsamples_=BitUtils::get32LH(buf);
  //cout << "numsamples: " << numsamples << endl;
  //cout << " filepos:   " << fin.file.tellg() << endl;
  std::vector <uint8_t>profile_buf(profile_size_bytes_);
  fin.file.read(reinterpret_cast<char*>(&profile_buf[0]),profile_size_bytes_);
  DecodeProfile(baseprofile,profile_buf);

  for (int ch=0;ch<numchannels_;ch++) {
    fin.file.read(reinterpret_cast<char*>(buf),6);
    uint32_t blocksize=BitUtils::get32LH(buf);
    uint16_t flag=BitUtils::get16LH(buf+4);
    if (flag>>9) framestats[ch].enc_mapped=true;
    else framestats[ch].enc_mapped=false;
    framestats[ch].maxbpn=flag&0xff;
    fin.ReadData(encoded[ch].GetBuf(),blocksize);
  }
}

void Codec::PrintProgress(int samplesprocessed,int totalsamples)
{
  double r=samplesprocessed*100.0/(double)totalsamples;
  std::cout << "  " << samplesprocessed << "/" << totalsamples << ":" << std::setw(6) << miscUtils::ConvertFixed(r,1) << "%\r";
}

void Codec::ScanFrames(Sac &mySac)
{
  std::streampos fsize=mySac.getFileSize();
  int frame_num=1;
  while (mySac.file.tellg()<fsize) {
    uint8_t buf[8];
    mySac.file.read(reinterpret_cast<char*>(buf),4);
    int numsamples=BitUtils::get32LH(buf);
    std::cout << "Frame " << frame_num << ": " << numsamples << " samples "<< std::endl;

    for (int ch=0;ch<mySac.getNumChannels();ch++) {
      mySac.file.read(reinterpret_cast<char*>(buf),6);
      uint32_t blocksize=BitUtils::get32LH(buf);
      uint16_t flag=BitUtils::get16LH(buf+4);
      //bool enc_mapped=flag>>9;
      int maxbpn=flag&0xff;
      std::cout << "  Channel " << ch << ": " << blocksize << " bytes (bpn: " << maxbpn << ", flags: " << (flag>>9) << ")" << std::endl;
      mySac.file.seekg(blocksize,std::ios_base::cur);
    }
    frame_num++;
  }
}

void Codec::EncodeFile(Wav &myWav,Sac &mySac,int profile)
{
  const int numchannels=myWav.getNumChannels();
  framesize=4*myWav.getSampleRate();

  FrameCoder myFrame(numchannels,framesize,profile,0);

  mySac.SetProfile(profile);
  mySac.WriteHeader(myWav);
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
  gtimer.stop();
  double time_total=gtimer.elapsedS();
  if (time_total>0.)   {
     double rprd=time_prd*100./time_total;
     double renc=time_enc*100./time_total;
     std::cout << "Timing: pred " << miscUtils::ConvertFixed(rprd,2) << "%, ";
     std::cout << "enc " << miscUtils::ConvertFixed(renc,2) << "%, ";
     std::cout << "misc " << miscUtils::ConvertFixed(100.-rprd-renc,2) << "%" << std::endl;
  }
}

void Codec::DecodeFile(Sac &mySac,Wav &myWav)
{
  const int numchannels=mySac.getNumChannels();
  framesize=4*mySac.getSampleRate();

  int samplestodecode=mySac.getNumSamples();
  int samplesdecoded=0;

  myWav.InitFileBuf(framesize);
  mySac.UnpackMetaData(myWav);
  myWav.WriteHeader();
  int profile=mySac.GetProfile();
  FrameCoder myFrame(numchannels,framesize,profile,0);
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
  std::cout << std::endl;
}
