#ifndef DE_H
#define DE_H

#include "opt.h"

#define DE_GENF_JADE

// Differential Evolution
class OptDE : public Opt {
  public:
    enum MutMethod {BEST1BIN,RAND1BIN,CUR1BEST};
    enum InitMethod {INIT_UNIV,INIT_NORM};
    struct DECfg
    {
      int NP=30;
      double CR=0.5;
      double F=0.5;
      double c=0.1;
      MutMethod mut_method=CUR1BEST;
      InitMethod init_method=INIT_NORM;
      double sigma_init=0.15;
      std::size_t num_threads=1;
      std::size_t nfunc_max=0;
    };

    OptDE(const DECfg &cfg,const box_const &parambox,bool verbose=false);

    ppoint run(opt_func func,const vec1D &xstart) override;


  protected:
    auto generate_candidate(const opt_points &pop,const vec1D &xbest,int iagent,double mCR,double mF);

    double gen_CR(double mCR)
    {
      return std::clamp(rand.r_norm(mCR,0.1),0.01,1.0);
    }
    double gen_F(double mF)
    {
      #ifdef DE_GENF_JADE
        double t=rand.event(1./3.)?
          rand.r_int(0,1.2):rand.r_norm(mF,0.1);
      #else
        double t=rand.r_norm(mF,0.1);
      #endif

      return std::clamp(t,0.01,1.2);
    }

    vec1D mut_1bin(const vec1D &xb,const vec1D &x1,const vec1D &x2,double F);
    vec1D mut_curbest(const vec1D &xbest,const vec1D &xb,const vec1D &x1,const vec1D &x2,double F);

    const DECfg &cfg;
    bool verbose;
};

#endif // DDS_H


