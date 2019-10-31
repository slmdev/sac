#ifndef LMS_CASCADE_H
#define LMS_CASCADE_H

#include "lms.h"
#include "lpc.h"
#include "rls.h"

class LMSCascade {
  public:
    LMSCascade(const std::vector<int> &vn,const std::vector<double>&vmu,const std::vector<double>vmudecay,const std::vector<double> &vpowdecay,double mu_mix,double mu_mix_beta)
    :n(vn.size()),p(n),clms(n),lad_mix(n,mu_mix,mu_mix_beta)
    {
      for (int i=0;i<n;i++)
        clms[i]=new NLMS_ROLL(vn[i],vmu[i],vmudecay[i],vpowdecay[i]);
    }
    double Predict()
    {
      for (int i=0;i<n;i++) p[i]=clms[i]->Predict();
      return lad_mix.Predict(p);
    }
    void Update(double val)
    {
      double pr=0.0;
      for (int i=0;i<n; i++) {
        clms[i]->Update(val-pr);
        pr+=p[i];
      }
      lad_mix.Update(val);

    }
    ~LMSCascade()
    {
      for (int i=0;i<n;i++) delete clms[i];
    }
  private:
    int n;
    vec1D p;
    std::vector<NLMS_ROLL*> clms;
    LAD_ADA lad_mix;
};

#endif // LMS_CASCADE_H
