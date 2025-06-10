#ifndef LMS_CASCADE_H
#define LMS_CASCADE_H

#include "lms.h"
#include "rls.h"
#include "../common/utils.h"

// increase stability
constexpr int LMS_MIX_INIT=1;
constexpr int LMS_MIX_CLAMPW=1;

class LMSCascade {
  public:
    LMSCascade(const std::vector<int> &vn,const std::vector<double>&vmu,
               const std::vector<double>&vmudecay,const std::vector<double> &vpowdecay,
               double mu_mix,double mu_mix_beta,int lm_n,double lm_alpha)
    :n(vn.size()),p(n+1),
     mix(n+1,mu_mix,mu_mix_beta),
     lm(lm_n,lm_alpha),
     clms(n)
    {
      if constexpr(LMS_MIX_INIT)
        for (int i=0;i<n;i++) mix.w[i] = 1.0/(i+1);

      for (int i=0;i<n;i++)
        clms[i]=new NLMS_Stream(vn[i],vmu[i],vmudecay[i],vpowdecay[i]);
    }
    double Predict()
    {
      for (int i=0;i<n;i++)
          p[i]=clms[i]->Predict();
      p[n]=lm.Predict();
      return mix.Predict(p);
    }
    void Update(const double target)
    {
      mix.Update(target);
      if constexpr(LMS_MIX_CLAMPW) {
        for (int i=0;i<=n;i++)
          mix.w[i]=std::max(mix.w[i],0.0);
      }

      double t=target;
      for (int i=0;i<n; i++) {
        clms[i]->Update(t);
        t-=mix.w[i]*p[i];
      }
      lm.UpdateHist(t);
    }
    ~LMSCascade()
    {
      for (int i=0;i<n;i++) delete clms[i];
    }
  private:
    int n;
    vec1D p;
    LAD_ADA mix;
    RLS lm;
    std::vector<LS_Stream*> clms;
};

#endif
