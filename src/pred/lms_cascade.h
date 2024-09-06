#ifndef LMS_CASCADE_H
#define LMS_CASCADE_H

#include "lms.h"
#include "../common/utils.h"

#define LMS_DECAY_EXP
//#define LMS_ADA

class LMSCascade {
  public:
    LMSCascade(const std::vector<int> &vn,const std::vector<double>&vmu,const std::vector<double>vmudecay,const std::vector<double> &vpowdecay,double mu_mix,double mu_mix_beta,double mix_nu)
    :n(vn.size()),p(n),clms(n),
    lms_mix(n,mu_mix,mu_mix_beta),pnu(n)
    {
      for (int i=0;i<n;i++) {
        #ifdef LMS_DECAY_EXP
          pnu[i]=exp(-(1.0-mix_nu)*i);
        #else
          pnu[i]=mix_nu;
        #endif
      }
      #ifdef LMS_ADA
        for (int i=0;i<n-1;i++)
          clms[i]=new NLMS_Stream(vn[i],vmu[i],vmudecay[i],vpowdecay[i]);
        clms[n-1]=new LMSADA_Stream(vn[n-1],vmu[n-1],0.97);
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
    void Update(double input, double ols_pred)
    {
      double target=input-ols_pred;
      lms_mix.Update(target);

      for (int i=0;i<n; i++) {
        clms[i]->Update(target);
        target-=pnu[i]*p[i];
      }
    }
    ~LMSCascade()
    {
      for (int i=0;i<n;i++) delete clms[i];
    }
  private:
    int n;
    vec1D p;
    std::vector<LS_Stream*> clms;
    LAD_ADA lms_mix;
    vec1D pnu;
};

#endif // LMS_CASCADE_H
