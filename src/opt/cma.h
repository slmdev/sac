#ifndef CMA_H
#define CMA_H

#include "opt.h"

class OptCMA : public Opt {
  public:
    struct CMACfg
    {
      int num_threads=1;
      int nfunc_max=0;
    };
    struct mparams
    {
      double d,p_target_succ,cp;
      double cc,ccov,ccovm,pthres;
      double psucc,sigma;
    } p;
    OptCMA(const CMACfg &cfg,const box_const &parambox,bool verbose=false);
    ppoint run(opt_func func,const vec1D &xstart) override;
  protected:
    const CMACfg &cfg;
    bool verbose;
};

#endif

