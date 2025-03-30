#ifndef DDS_H
#define DDS_H

#include "opt.h"
#include <future>
#include <cassert>

#define DDS_SIGMA_ADAPT

// Dynamical dimensioned search algorithm for computationally efficient watershed model calibration
// Tolson, Shoemaker 2007
class OptDDS : public Opt {
  public:
    struct DDSCfg
    {
      double sigma_start=0.; // gets overwritten
      double sigma_max=0.5;
      double sigma_min=0.05;
      int c_succ_max=3;
      int c_fail_max=50;
      int num_threads=1;
      int nfunc_max=0;
    };

    OptDDS(const DDSCfg &cfg,const box_const &parambox,bool verbose=false)
    :Opt(parambox),cfg(cfg),verbose(verbose)
    {
    }
    vec1D generate_candidate(const vec1D &x,int nfunc,int nfunc_max,double sigma)
    {
      std::vector <int>J; // select J of D variables
      double p=1.0-log(nfunc)/log(nfunc_max);

      for (int i=0;i<ndim;i++) {
        if (rand.event(p)) J.push_back(i);
      }
      // set empty? select random element
      if (!J.size()) J.push_back(rand.ru_int(0,ndim-1));

      // perturb decision variables
      vec1D xtest=x;
      for (auto k:J) {
        xtest[k]=gen_norm(xtest[k],pb[k],sigma);
        assert(xtest[k]>=pb[k].xmin && xtest[k]<=pb[k].xmax);
      }
      return xtest;
    }

    // sequential single threaded
    ppoint run_single(opt_func func,const vec1D &xstart)
    {
      int nfunc=1;
      ppoint xb{func(xstart),xstart};
      if (verbose) std::cout << xb.first << '\n';

      double sigma=cfg.sigma_start;
      #ifdef DDS_SIGMA_ADAPT
        int c_succ=0,c_fail=0;
      #endif

      while (nfunc<cfg.nfunc_max) {
        ppoint x_gen;
        x_gen.second=generate_candidate(xb.second,nfunc,cfg.nfunc_max,sigma);
        x_gen.first=func(x_gen.second);

        #ifndef DDS_SIGMA_ADAPT
          if (x_gen.first<xb.first)
            xb = x_gen;
        #else
          if (x_gen.first<xb.first) {
            xb=x_gen;

            c_succ+=1;
            c_fail=0;
          } else {
            c_fail+=1;
            c_succ=0;
          };

          if (c_succ >= cfg.c_succ_max) {
            sigma=std::min(2.0*sigma,cfg.sigma_max);
            c_succ=0;
          } else if (c_fail >= cfg.c_fail_max) {
            sigma=std::max(sigma/2.0,cfg.sigma_min);
            c_fail=0;
          }
          #endif
        nfunc++;
        if (verbose) std::cout << " DDS " << std::format("{:5}",nfunc) << ": " << std::format("{:0.4f}",xb.first) << " s=" << sigma << "\r";
      }
      return xb;
    }

    // multi-threaded variant
    ppoint run_mt(opt_func func,const vec1D &xstart)
    {

      ppoint xb{func(xstart),xstart}; // eval at initial solution

      if (verbose) std::cout << xb.first << '\n';

      double sigma=cfg.sigma_start;

      int nfunc=1;
      while (nfunc<cfg.nfunc_max) {
        int nthreads = std::min(cfg.nfunc_max-nfunc,cfg.num_threads);

        // generate candidates around current xbest
        opt_points x_gen(nthreads);
        for (int i=0;i<nthreads;i++) {
          x_gen[i].second=generate_candidate(xb.second,nfunc,cfg.nfunc_max,sigma);
          nfunc++;
        };
        //nfunc+=nthreads;

        eval_points_mt(func,std::span(x_gen));

        // select
        for (int i=0;i<nthreads;i++)
          if (x_gen[i].first<xb.first)
            xb = x_gen[i];

        if (verbose) std::cout << " DDS mt=" << nthreads << ": " << std::format("{:5}",nfunc) << ": " << std::format("{:0.4f}",xb.first) << " s=" << std::format("{:0.3f}",sigma) << "\r";
      }
      return xb;
    }

    ppoint run(opt_func func,const vec1D &xstart)
    {
      assert(pb.size()==xstart.size());

      ppoint pbest;
      if (cfg.num_threads<=1) pbest=run_single(func,xstart);
      else pbest=run_mt(func,xstart);

      if (verbose) std::cout << '\n';
      return pbest;
    }
  protected:
    const DDSCfg &cfg;
    bool verbose;
};

#endif // DDS_H

