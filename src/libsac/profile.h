#ifndef PROFILE_H
#define PROFILE_H

#include <vector>
#include "../pred/lpc.h"
#include "../pred/nlms.h"
#include "../pred/rls.h"
#include "map.h"

class SacProfile {
  public:
    struct FrameStats {
      int maxbpn,maxbpn_map;
      bool enc_mapped;
      int32_t minval,maxval,mean;
      Remap mymap;
    };
    struct coef {
      float vmin,vmax,vdef;
    };
      SacProfile():type(0){};
      void Init(int numcoefs,int ptype)
      {
         coefs.resize(numcoefs);
         type=ptype;
      }
      SacProfile(int numcoefs,int type)
      :type(type),coefs(numcoefs)
      {

      }
      void Set(int num,double vmin,double vmax,double vdef)
      {
        if (num>=0 && num< static_cast<int>(coefs.size())) {
          coefs[num].vmin=vmin;
          coefs[num].vmax=vmax;
          coefs[num].vdef=vdef;
        }
      }
      double Get(int num) const {
        if (num>=0 && num< static_cast<int>(coefs.size())) {
            return coefs[num].vdef;
        } else return 0.;
      }
      int type;
      std::vector <coef> coefs;
};

// abstract predictor class
class StereoPredictor {
  public:
    StereoPredictor():ch0(0),ch1(1){};
    virtual int32_t Predict()=0;
    virtual void Update(int32_t val)=0;
    virtual ~StereoPredictor(){};
  protected:
    int ch0,ch1;
};

class StereoFast : public StereoPredictor {
  public:
    StereoFast(const SacProfile &profile);
    int32_t Predict();
    void Update(int32_t val);
  private:
    double x[2][2],p[2];
    std::vector <LPC2> lpc;
    std::vector <NLMS> nlms;
};

class StereoNormal : public StereoPredictor {
  static const int nstages=4;
  public:
    StereoNormal(const SacProfile &profile);
    int32_t Predict();
    void Update(int32_t val);
  private:
    double x[2][nstages],p[nstages];
    std::vector <LPC3> lpc;
    std::vector <LMSADA2> lms1,lms2,lms3;
    std::vector <RLS> lmsmix;
    std::vector<double> pv;
};

class StereoHigh : public StereoPredictor {
  static const int nstages=5;
  public:
    StereoHigh(const SacProfile &profile);
    int32_t Predict();
    void Update(int32_t val);
  private:
    double x[2][nstages],p[nstages];
    std::vector <LPC3> lpc;
    std::vector <LMSADA2>lms1,lms2,lms3,lms4;
    std::vector <RLS> lmsmix;
    std::vector<double> pv;
};


#endif // PROFILE_H
