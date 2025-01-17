#ifndef DE_H
#define DE_H

#include "opt.h"
#include "dds.h"
#include <cassert>

#define DE_GENF_JADE

// Differential Evolution
class OptDE : public Opt {
  public:
    enum MutMethod {BEST1BIN,RAND1BIN,CUR1BEST};
    enum InitMethod {INIT_RAND,INIT_NORM};
    struct DECfg
    {
      int NP=50;
      double CR=0.5;
      double F=0.6;
      double c=0.1;
      MutMethod mut_method=CUR1BEST;
      InitMethod init_method=INIT_RAND;
      double init_sigma=0.15;
    };
    OptDE(const box_const &parambox)
    :Opt(parambox),ndim(pb.size())
    {
    }

    ppoint run(const DECfg &cfg,opt_func func,const vec1D &xstart,int nfunc_max,bool verbose=true)
    {
      assert(pb.size()==xstart.size());
      std::vector<ppoint> pop(cfg.NP);

      int nfunc=0;
      vec1D xbest(ndim);
      double fbest=std::numeric_limits<double>::max();

      if (cfg.init_method == INIT_RAND || cfg.init_method == INIT_NORM)
      {
        xbest=xstart;
        fbest=func(xstart); // eval at initial solution
        pop[0] = {fbest,xbest};
        nfunc++;

        for (int i=1;i<cfg.NP;i++) { // fill population
          vec1D x;
          if (cfg.init_method == INIT_RAND) x=gen_uniform_sample(xbest,cfg.init_sigma);
          else if (cfg.init_method == INIT_NORM)x=gen_norm_sample(xbest,cfg.init_sigma);

          double fx=func(x); // eval
          nfunc++;

          pop[i] = {fx,x};

          if (fx < fbest) {
            fbest = fx;
            xbest = x;
          }
        }
      } else std::cerr << "INIT_DDS: POP not initialized correctly\n";

      if (verbose) std::cout << "\n" << " NP=" << cfg.NP << " " << fbest << "\n";

      double mCR = cfg.CR;
      double mF = cfg.F;

      while (nfunc<nfunc_max)
      {
        std::vector<double>CR_succ;
        std::vector<double>F_succ;

        for (size_t iagent=0;iagent<pop.size();++iagent)
        {
          const double tCR = gen_CR(mCR);
          const double tF  = gen_F(mF);
          const int R = rand.ru_int(0,ndim-1);

          // mutation
          std::vector<int> e(pop.size());
          std::iota(std::begin(e),std::end(e),0);std::erase(e,iagent);

          vec1D xm;
          if (cfg.mut_method==BEST1BIN) {
            int xr1 = e[rand.ru_int(0,e.size()-1)];std::erase(e,xr1);
            int xr2 = e[rand.ru_int(0,e.size()-1)];std::erase(e,xr2);
            xm = mut_1bin(xbest,pop[xr1].second,pop[xr2].second,mF);
          } else if (cfg.mut_method==RAND1BIN) {
            int xr0 = e[rand.ru_int(0,e.size()-1)];std::erase(e,xr0);
            int xr1 = e[rand.ru_int(0,e.size()-1)];std::erase(e,xr1);
            int xr2 = e[rand.ru_int(0,e.size()-1)];std::erase(e,xr2);
            xm = mut_1bin(pop[xr0].second,pop[xr1].second,pop[xr2].second,mF);
          } else if (cfg.mut_method==CUR1BEST) {
            int xr1 = e[rand.ru_int(0,e.size()-1)];std::erase(e,xr1);
            int xr2 = e[rand.ru_int(0,e.size()-1)];std::erase(e,xr2);
            xm = mut_curbest(xbest,pop[iagent].second,pop[xr1].second,pop[xr2].second,mF);
          }

          // cross-over
          vec1D xtrial(ndim);
          const ppoint &xi=pop[iagent];
          for (int i=0;i<ndim;i++)   {
            if (rand.event(mCR) || (i==R)) xtrial[i] = xm[i];
            else xtrial[i] = xi.second[i];
          }
          double fx=func(xtrial);
          nfunc++;

          // greedy selection
          if (fx < xi.first) {
            pop[iagent] = {fx,xtrial}; // replace
            CR_succ.push_back(tCR);
            F_succ.push_back(tF);
          }
          if (verbose) std::cout << "DE " << nfunc << ": " << fbest << " (mCR=" << mCR << " mF=" << mF << ")\r";

          if (fx < fbest) {
            fbest = fx;
            xbest = xtrial;
          }

          if (nfunc >= nfunc_max) break;
        }
        mCR = (1.0-cfg.c)*mCR + cfg.c*MathUtils::mean(CR_succ);
        mF  = (1.0-cfg.c)*mF  + cfg.c*MathUtils::mean(F_succ);
      }
      if (verbose) std::cout << '\n';
      return {fbest,xbest};
    }
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

    vec1D mut_1bin(const vec1D &xb,const vec1D &x1,const vec1D &x2,double F)
    {
      vec1D xm(ndim);
      for (int i=0;i<ndim;i++) {
        double y=xb[i] + F*(x1[i] - x2[i]);
        xm[i] = reflect(y, pb[i].xmin,pb[i].xmax);
      }
      return xm;
    }
    vec1D mut_curbest(const vec1D &xbest,const vec1D &xb,const vec1D &x1,const vec1D &x2,double F)
    {
      vec1D xm(ndim);
      for (int i=0;i<ndim;i++) {
        double y=xb[i] + F*(xbest[i] - xb[i]) + F*(x1[i] - x2[i]);
        xm[i] = reflect(y, pb[i].xmin,pb[i].xmax);
      }
      return xm;
    }
    int ndim;
};

#endif // DDS_H


