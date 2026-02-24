#ifndef LMS_CASCADE_H
#define LMS_CASCADE_H

#include "lms.h"
#include "rls.h"
#include "blend.h"
#include "../common/utils.h"




static std::vector<std::unique_ptr<LMS>>
make_mix(int n,double mu_mix,double mu_mix_beta)
{
    std::vector<std::unique_ptr<LMS>> p;
    p.push_back(std::make_unique<LAD_ADA>(n, 1*mu_mix,mu_mix_beta));
    p.push_back(std::make_unique<LMS_ADA>(n, 1*mu_mix,mu_mix_beta));
    //p.push_back(std::make_unique<LAD_ADA>(n, 2*mu_mix,mu_mix_beta));
    //p.push_back(std::make_unique<LMS_ADA>(n, 2*mu_mix,mu_mix_beta));

    if constexpr(SACCfg::LMS_MIX_INIT) {
      for (std::size_t i=0;i<p.size();i++)
        for (int k=0;k<n;k++) p[i]->w[k] = 1.0/(1.0+k);
    }

    return p;
}

class Cascade {
  public:
    Cascade(const std::vector<int> &vn,const std::vector<double>&vmu,
               const std::vector<double>&vmudecay,const std::vector<double> &vpowdecay,
               double mu_mix,double mu_mix_beta,int lm_n,double lm_alpha)
    :n(vn.size()),p(n+1),
     mix(make_mix(n+1,mu_mix,mu_mix_beta)),
     lm(lm_n,lm_alpha),
     clms(n)
    {
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

      double t=target;
      for (int i=0;i<n; i++) {
        clms[i]->Update(t);
        t-=mix.GetWeightedContribution(i)*p[i];
      }
      lm.UpdateHist(t);

      mix.Update(target);
    }
    ~Cascade()
    {
      for (int i=0;i<n;i++) delete clms[i];
    }
  private:
    int n;
    vec1D p;
    BlendLMS mix;
    RLS lm;
    std::vector<LS_Stream*> clms;
};

#endif
