#ifndef CASCADE_H
#define CASCADE_H

#include <memory>
#include "blend.h"
#include "ls.h"
#include "rls.h"
#include "../common/utils.h"

// horizontal blend of arbitrary number of (LS-)experts
class BlendLS {
  public:
    explicit BlendLS(std::vector<std::unique_ptr<LS>> preds,double mv_alpha=0.9,double mv_beta=1.0)
    :expert(std::move(preds)), //pp has ownership of predictors
    np(expert.size()),
    ep(np),ew(np),
    sm(np,mv_alpha,mv_beta)
    {
      if (expert.empty())
            throw std::invalid_argument("BlendLS: preds.size() must be >0");
    }

    //blended contribution of index component over all experts
    double GetWeight(int index) {
        for (std::size_t i=0;i<np;i++)
          ew[i]=expert[i]->GetWeight(index);

        double weight=slmath::dot(ew,sm.Weights());
        return std::max(weight,0.0);
    }
    //prediction using blended predictors over input
    double Predict(const vec1D &input)
    {
      for (std::size_t i=0;i<np;i++) {
        ep[i]=expert[i]->Predict(input); // call expert i
        if (!std::isfinite(ep[i]))
          throw std::invalid_argument("exception in cascade, predictor is not finite");
      }
      return sm.Predict(ep); //blend;
    }
    void Update(double target)
    {
      for (std::size_t i=0;i<np;i++)
        expert[i]->Update(target);
      sm.Update(target);
    }

private:
    std::vector<std::unique_ptr<LS>> expert;
    std::size_t np;
    vec1D ep,ew;
    BlendExp<RunSumEMA> sm;
};


static std::vector<std::unique_ptr<LS>>
make_mix(int n,double mu_mix,double mu_mix_beta)
{
    std::vector<std::unique_ptr<LS>> p;
    p.push_back(std::make_unique<LS_ADA<Loss::L1,LSInitType::Uniform>>(n,mu_mix,mu_mix_beta));
    p.push_back(std::make_unique<LS_ADA<Loss::L2,LSInitType::Uniform>>(n,mu_mix,mu_mix_beta));
    //p.push_back(std::make_unique<HM::HMix<Loss::L1,HM::Gate1<HM::Tanh>,HM::NoReg,2>>(n,mu_mix,mu_mix_beta));
    //p.push_back(std::make_unique<HM::HMix<Loss::L2,HM::Gate1<HM::Tanh>,HM::NoReg,2>>(n,mu_mix,mu_mix_beta));
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
        t-=mix.GetWeight(i)*p[i];
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
    BlendLS mix;

    RLS lm;
    std::vector<LS_Stream*> clms;
};

#endif
