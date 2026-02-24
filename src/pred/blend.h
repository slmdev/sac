#ifndef BLEND_H
#define BLEND_H

#include <memory>
#include "lms.h"

// blend two expert outputs via sigmoid
class Blend2 {
  public:
    Blend2(double beta=0.95,double theta0=1.0,double theta1=0.0,double scale=5.0)
    :w(0.5),
     th0(theta0),th1(theta1),scale(scale),
     rsum(beta)
    {
    }
    double Predict(double p0,double p1) const
    {
      return w*p0 + (1.0-w)*p1;
    }

    void Update(double score0,double score1)
    {
      rsum.Update(score1-score0);

      double delta=rsum.Get();
      double z = th0*delta + th1;
      z = std::clamp(z,-scale,scale);
      w = 1.0 / (1.0 + std::exp(-z));
   }
  private:
    double w,th0,th1,scale;
    RunSumEMA rsum;
};

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
    pv(np),zm(np),w(np),
    rsum(np,RunMeanVar(mv_alpha)),beta(mv_beta),px(0.0)
    {
      if (expert.empty())
            throw std::invalid_argument("BlendLMS: preds.size() must be >0");
      std::fill(begin(w),end(w),1.0/np); //init equal weight
    }

    //blended contribution of index component over all experts
    double GetWeightedContribution(int index) const {

        double weight=0.0;
        for (std::size_t i=0;i<np;i++) {
          weight=std::fma(w[i],expert[i]->w[index],weight);
        }
        return std::max(weight,0.0);
    }
    //prediction using blended predictors over input
    double Predict(const vec1D &input)
    {
      for (std::size_t i=0;i<np;i++) {
        pv[i]=expert[i]->Predict(input); // call expert i
        if (!std::isfinite(pv[i]))
          throw std::invalid_argument("exception in cascade, predictor is not finite");
      }

      px=slmath::dot(w,pv); //blend
      return px;
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

      //update blending weights using old prediction
      double loss_px = std::abs(target-px);
      for (std::size_t i=0;i<np;i++) {
        double loss_pi=std::abs(target-pv[i]);
        // if regret < 0 -> expert better then blend
        double regret=(loss_pi - loss_px);
        rsum[i].Update(regret);
      }

      #if 0
        double hr=0.5,beta_nu=0.05;
        const double Hmax = std::log(np);
        const double Htarget = Hmax*hr;
        double H = 0.0;
        for (std::size_t i=0;i<np;i++) {
          if (w[i]>1E-8) H-=w[i]*log(w[i]);
        }
        beta += beta_nu*(H-Htarget);
        beta = std::clamp(beta,1.0,15.0);
        //std::cout << H << ' ' << Htarget << ' ' << beta << '\n';
      #endif

      UpdateWeights(w,beta);
    }

private:
    // softmax w_i = exp(-beta * normalized_regret)
    void UpdateWeights(vec1D &wd,double beta)
    {
      constexpr double EPS=1E-8;
      double max_z = -std::numeric_limits<double>::infinity();
      for (std::size_t i=0;i<np;i++) {
        auto [mean,var] = rsum[i].get(); //regret

        // scaled signal-to-noise
        zm[i]= -beta*mean/(std::sqrt(var)+EPS);
        max_z = std::max(max_z,zm[i]);
      }

      //best expert has highest z-score -> weight=exp(0)=1
      double total=0.0;
      for (std::size_t i=0;i<np;i++) {
        wd[i] = std::exp(zm[i]-max_z);
        total += wd[i];
      }
      //normalize weights, total >= 1 from max-trick
      const double inv_total=1.0/total;
      for (double &val : wd) val *= inv_total;
    }

    std::vector<std::unique_ptr<LMS>> expert;
    std::size_t np;
    vec1D pv,zm,w;
    std::vector <RunMeanVar> rsum;
    double beta;
    double px;
};


#endif

