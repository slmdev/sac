#include "profile.h"

StereoFast::StereoFast(const SacProfile &profile)
:lpc(2,LPC2(profile.Get(0),16,8,32)),
nlms(2,NLMS(profile.Get(1),64))
{
  for (int i=0;i<2;i++) x[0][i]=x[1][i]=0.;
};

int32_t StereoFast::Predict()
{
  p[0]=lpc[ch0].Predict(x[ch1][0]);
  p[1]=nlms[ch0].Predict();
  return std::round(p[0]+p[1]);
}

void StereoFast::Update(int32_t val)
{
  double *px=x[ch0];

  px[0]=val;
  px[1]=px[0]-p[0];

  lpc[ch0].Update(px[0]);
  nlms[ch0].Update(px[1]);
  std::swap(ch0,ch1);
}

StereoNormal::StereoNormal(const SacProfile &profile)
:lpc(2,LPC3(profile.Get(0),16,8,32,500)),
 lms1(2,LMSADA2(256,128,profile.Get(1),profile.Get(2))),
 lms2(2,LMSADA2(32,16,profile.Get(3),profile.Get(3))),
 lms3(2,LMSADA2(8,4,profile.Get(4),profile.Get(4))),
 lmsmix(2,RLS(3,profile.Get(5))),
 pv(3)
{
  for (int i=0;i<nstages;i++) x[0][i]=x[1][i]=0;
};

int32_t StereoNormal::Predict()
{
  double *px=x[ch1];

  p[0]=lpc[ch0].Predict(px[0]);
  p[1]=lms1[ch0].Predict(px[1]);
  p[2]=lms2[ch0].Predict(px[2]);
  p[3]=lms3[ch0].Predict(px[3]);

  pv[0]=p[1];pv[1]=pv[0]+p[2];pv[2]=pv[1]+p[3];
  return round(p[0]+lmsmix[ch0].Predict(pv));
}

void StereoNormal::Update(int32_t val)
{
  double *px=x[ch0];

  px[0]=val;
  px[1]=px[0]-p[0];
  px[2]=px[1]-p[1];
  px[3]=px[2]-p[2];

  lpc[ch0].Update(px[0]);
  lms1[ch0].Update(px[1]);
  lms2[ch0].Update(px[2]);
  lms3[ch0].Update(px[3]);
  lmsmix[ch0].Update(px[1]);
  std::swap(ch0,ch1);
}

StereoHigh::StereoHigh(const SacProfile &profile)
:lpc(2,LPC3(profile.Get(0),16,8,8,500)),
 lms1(2,LMSADA2(1280,640,profile.Get(1),profile.Get(2))),
 lms2(2,LMSADA2(256,128,profile.Get(3),profile.Get(4))),
 lms3(2,LMSADA2(32,16,profile.Get(5),profile.Get(5))),
 lms4(2,LMSADA2(8,4,profile.Get(6),profile.Get(6))),
 lmsmix(2,RLS(4,profile.Get(7))),
     //lmsmix(2,BlendRPROP(4,0.005,0.001,0.01)),
 pv(4)
{
  for (int i=0;i<5;i++) x[0][i]=x[1][i]=0;
};

int32_t StereoHigh::Predict()
{
  p[0]=lpc[ch0].Predict(x[ch1][0]);
  p[1]=lms1[ch0].Predict(x[ch1][1]);
  p[2]=lms2[ch0].Predict(x[ch1][2]);
  p[3]=lms3[ch0].Predict(x[ch1][3]);
  p[4]=lms4[ch0].Predict(x[ch1][4]);

  pv[0]=p[1];pv[1]=pv[0]+p[2];pv[2]=pv[1]+p[3];pv[3]=pv[2]+p[4];
  return round(p[0]+lmsmix[ch0].Predict(pv));
}

void StereoHigh::Update(int32_t val)
{
  x[ch0][0]=val;
  x[ch0][1]=x[ch0][0]-p[0];
  x[ch0][2]=x[ch0][1]-p[1];
  x[ch0][3]=x[ch0][2]-p[2];
  x[ch0][4]=x[ch0][3]-p[3];

  lpc[ch0].Update(x[ch0][0]);
  lms1[ch0].Update(x[ch0][1]);
  lms2[ch0].Update(x[ch0][2]);
  lms3[ch0].Update(x[ch0][3]);
  lms4[ch0].Update(x[ch0][4]);
  lmsmix[ch0].Update(x[ch0][1]);
  std::swap(ch0,ch1);
}

