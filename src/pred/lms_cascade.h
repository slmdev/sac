#ifndef LMS_CASCADE_H
#define LMS_CASCADE_H

#include <memory>
#include "blend.h"
#include "lms.h"
#include "rls.h"
#include "../common/utils.h"

// Blend 2xLMS-ADA using L1 + L2 loss
// using absolute error as scoring function
class Blend2LMS_L1 {
  public:
    Blend2LMS_L1(int n,double lms_mu,double lms_beta)
    :n(n),px0(0.0),px1(0.0),
     mix0(n,lms_mu,lms_beta),
     mix1(n,lms_mu,lms_beta),
     cw2(lms_beta)
    {
      if constexpr(SACCfg::LMS_MIX_INIT)
        for (int i=0;i<n-1;i++)
          mix0.w[i] = mix1.w[i] = 1.0/(i+1);
    }
    double GetWeight(int index) const
    {
      return cw2.Predict(mix0.w[index],mix1.w[index]);
    }
    double Predict(const vec1D &input)
    {
      px0=mix0.Predict(input);
      px1=mix1.Predict(input);
      return cw2.Predict(px0,px1);
    }
    void UpdateMixer(double target)
    {
      mix0.Update(target);
      mix1.Update(target);
      if constexpr(SACCfg::LMS_MIX_CLAMPW)
        for (int i=0;i<n;i++) {
          mix0.w[i]=std::max(mix0.w[i],0.0);
          mix1.w[i]=std::max(mix1.w[i],0.0);
        }
    }
    void UpdateBlend(double target)
    {
      double e0=std::abs(target-px0);
      double e1=std::abs(target-px1);
      cw2.Update(e0,e1);
    }
    int n;
    double px0,px1;
    LAD_ADA mix0;
    LMS_ADA mix1;
    Blend2 cw2;
};


// horizontal blend of arbitrary number of (LMS-)experts
class BlendLMS {
  public:
    explicit BlendLMS(std::vector<std::unique_ptr<LMS>> preds,double mv_alpha=0.9,double mv_beta=5)
    :expert(std::move(preds)), //pp has ownership of predictors
    np(expert.size()),
    ep(np),ew(np),
    //rsum(np,RunMeanVar(mv_alpha)),beta(mv_beta),
    sm(np,mv_alpha,mv_beta)
    {
      if (expert.empty())
            throw std::invalid_argument("BlendLMS: preds.size() must be >0");
    }

    //blended contribution of index component over all experts
    double GetWeightedContribution(int index) {
        for (std::size_t i=0;i<np;i++)
          ew[i]=expert[i]->w[index];

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
      for (std::size_t i=0;i<np;i++) {
        expert[i]->Update(target);

        if constexpr(SACCfg::LMS_MIX_CLAMPW) {
          for (int k=0;k<(int)expert[i]->w.size();k++)
            expert[i]->w[k]=std::max(expert[i]->w[k],0.0);
        }
      }

      sm.Update(target);
    }

private:
    std::vector<std::unique_ptr<LMS>> expert;
    std::size_t np;
    vec1D ep,ew;
    //std::vector <RunMeanVar> rsum;
    //double beta;
    //double px;
    BlendExp<RunMeanVar> sm;
};


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
