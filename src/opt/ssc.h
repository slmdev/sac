#ifndef SSC_H
#define SSC_H

// adaptive step-size control

class SSC0
{
  public:
    SSC0(int nsucc_max=3,int nfail_max=50)
    :nsucc_max(nsucc_max),nfail_max(nfail_max)
    {
      nsucc = nfail=0;
    }
    double update(double sigma,double lambda)
    {
      if (lambda > 0.0) {
        nsucc+=1;
        nfail=0;
      } else {
        nsucc=0;
        nfail+=1;
      }

      if (nsucc >= nsucc_max) {
        sigma=sigma*2.0;
        nsucc=0;
      } else if (nfail >= nfail_max) {
        sigma=sigma/2.0;
        nfail=0;
      }
      return std::clamp(sigma,0.05,0.5);
    }
  protected:
    int nsucc_max,nfail_max,nsucc,nfail;
};

// p_target: target succession prob.
// p_a: learning rate for exp smoothed succession prob.
// r_sigma: relative increase/decrease for sigma
class SSC1
{
  public:
    SSC1(double p_target=0.10,double p_c=0.10,double p_d=0.05,double lambda_min=0.05,double lambda_max=0.25)
    :p_target(p_target),p_c(p_c),p_d(p_d),
    lmin(lambda_min),lmax(lambda_max)
    {
      p_succ = p_target;
    }
    double update(double sigma,double lambda)
    {
      // update empirical success prob by exp. smoothing
      p_succ=(1.0-p_c)*p_succ + p_c*lambda;
      sigma = sigma * std::exp(p_d * (p_succ-p_target) / (1.0-p_target));
      return std::clamp(sigma,lmin,lmax);
    }
    double p_succ;
  protected:
    double p_target,p_c,p_d;
    double lmin,lmax;
};

#if 0
  /*if (p_succ > p_target)
    sigma = std::min(sigma*1.01,0.5);
  else if (p_succ < p_target) {
    sigma = std::max(sigma/1.01,0.05);
  }*/
  double pd=0.01;
  sigma = sigma * std::exp(pd * (p_succ-p_target) / (1.0-p_target));
  sigma = std::max(std::min(sigma,0.5),0.05);
#endif


#endif // header guard

