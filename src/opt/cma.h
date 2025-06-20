#ifndef CMA_H
#define CMA_H

#include "opt.h"
#include "../common/utils.h"
#include "../common/math.h"

class OptCMA : public Opt {
  public:
    struct CMACfg
    {
      int num_threads=1;
      int nfunc_max=0;
      double sigma_init=0.;
    };
    struct CMAParams
    {
      double d,p_target_succ,cp;
      double cc,ccov,ccovm,pthres;
      double psucc,sigma;

      CMAParams(int n)
      {
        d=1.0+n/2.0; // damping parameter
        p_target_succ=2.0/11.0;
        cp=1.0/12.0; // learning rate

        // covariance matrix adaption
        cc=2.0/(n+2.0); // learning rate evolution path
        ccov=2.0/(n*n + 6.0); // learning rate covariance matrix
        ccovm=0.4/(std::pow(n,1.6)+1.0); // learning rate active covariance update
        pthres=0.44;

        psucc=sigma=0.0;
      }
    };
    OptCMA(const CMACfg &cfg,const box_const &parambox,bool verbose=false);
    ppoint run(opt_func func,const vec1D &xstart) override;
  protected:
    auto generate_candidate(const vec1D &x,double sigma);
    void update_cov(vec2D &mcov, vec1D &pc,const vec1D &az);
    const CMACfg &cfg;
    CMAParams p;
    slmath::Cholesky chol;
    bool verbose;

};

#endif

