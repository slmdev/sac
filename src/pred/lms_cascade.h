#ifndef LMS_CASCADE_H
#define LMS_CASCADE_H

#include "lms.h"
#include "../common/utils.h"

//#define LMS_ADA
//#define LMS_N0
#define LMS_INIT

class LMSCascade {
  public:
    LMSCascade(const std::vector<int> &vn,const std::vector<double>&vmu,const std::vector<double>vmudecay,const std::vector<double> &vpowdecay,double mu_mix,double mu_mix_beta)
    :n(vn.size()),
    #ifdef LMS_N0
      p(n+1),lms_mix(n+1,mu_mix,mu_mix_beta),
    #else
      p(n),lms_mix(n,mu_mix,mu_mix_beta),
    #endif
    clms(n)
    {
      #ifdef LMS_INIT
        for (int i=0;i<n;i++) lms_mix.w[i] = 1.0/(i+1);
      #endif
      #ifdef LMS_ADA
        for (int i=0;i<n-1;i++)
          clms[i]=new NLMS_Stream(vn[i],vmu[i],vmudecay[i],vpowdecay[i]);
        clms[n-1]=new LMSADA_Stream(vn[n-1],vmu[n-1],vmudecay[n-1],vpowdecay[n-1]);
      #else
        for (int i=0;i<n;i++)
          clms[i]=new NLMS_Stream(vn[i],vmu[i],vmudecay[i],vpowdecay[i]);
      #endif
    }
    double Predict()
    {
      for (int i=0;i<n;i++) {
          p[i]=clms[i]->Predict();
      }
      return lms_mix.Predict(p);
    }
    void Update(const double target)
    {
      lms_mix.Update(target);

      double t=target;
      for (int i=0;i<n; i++) {
        clms[i]->Update(t);
        t-=lms_mix.w[i]*p[i];
      }
      #ifdef LMS_N0
        p[n] = t;
      #endif
    }
    ~LMSCascade()
    {
      for (int i=0;i<n;i++) delete clms[i];
    }
  private:
    int n;
    vec1D p;
    LAD_ADA lms_mix;
    std::vector<LS_Stream*> clms;
};

#endif // LMS_CASCADE_H
