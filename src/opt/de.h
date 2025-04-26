#ifndef DE_H
#define DE_H

#include "opt.h"

// Differential Evolution
class OptDE : public Opt {
  public:
    enum MutMethod {BEST1BIN,RAND1BIN,CUR1BEST,CURPBEST};
    std::unordered_map<MutMethod, int> MutVals = {
      {BEST1BIN,2},
      {RAND1BIN,3},
      {CUR1BEST,2},
      {CURPBEST,2}
    };
    enum InitMethod {INIT_UNIV,INIT_NORM};
    struct DECfg
    {
      int NP=30;
      double CR=0.5;
      double F=0.5;
      double c=0.1;
      MutMethod mut_method=CURPBEST;
      InitMethod init_method=INIT_NORM;
      double pbest=0.1;
      int npbest=std::clamp(static_cast<int>(std::round(pbest*NP))-1,0,NP-1);
      double sigma_init=0.15;
      std::size_t num_threads=1;
      std::size_t nfunc_max=0;
    };

    OptDE(const DECfg &cfg,const box_const &parambox,bool verbose=false);

    ppoint run(opt_func func,const vec1D &xstart) override;


  protected:
    auto generate_candidate(const opt_points &pop,const vec1D &xbest,int iagent,double mCR,double mF);
    void print_status(std::size_t nfunc,double fx,double mCR,double mF);
    std::vector<int> select_k_unique_except(int n,int t,int k);

    double gen_CR(double mCR)
    {
      return std::clamp(rand.r_norm(mCR,0.1),0.01,1.0);
    }
    double gen_F(double mF)
    {
      return std::clamp(rand.r_cauchy(mF,0.1),0.01,1.0);
    }

    vec1D mut_1bin(const vec1D &xb,const vec1D &x1,const vec1D &x2,double F);
    vec1D mut_curbest(const vec1D &xbest,const vec1D &xb,const vec1D &x1,const vec1D &x2,double F);

    const DECfg &cfg;
    bool verbose;
};

#endif // DDS_H


