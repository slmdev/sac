#ifndef DDS_H
#define DDS_H

#include "opt.h"
#include <cassert>

#define DDS_SIGMA_ADAPT

// Dynamical dimensioned search algorithm for computationally efficient watershed model calibration
// Tolson, Shoemaker 2007
class OptDDS : public Opt {
  public:
    struct DDSCfg
    {
      double sigma_start=0.2;
      double sigma_max=0.5;
      double sigma_min=0.05;
      int c_succ_max=3;
      int c_fail_max=50;
    };

    OptDDS(const box_const &parambox)
    :Opt(parambox),ndim(parambox.size())
    {
    }
    ppoint run(const DDSCfg &cfg,opt_func func,const vec1D &xstart,int nfunc_max,bool verbose=true)
    {
      assert(pb.size()==xstart.size());

      int nfunc=1;
      vec1D xbest=xstart;
      double fbest=func(xstart); // eval at initial solution

      double sigma=cfg.sigma_start;

      #ifdef DDS_SIGMA_ADAPT
        int c_succ=0,c_fail=0;
      #endif

      if (verbose) std::cout << fbest << '\n';

      while (nfunc<nfunc_max) {
        std::vector <int>J; // select J of D variables
        double p=1.0-log(nfunc)/log(nfunc_max);

        for (int i=0;i<ndim;i++) {
          if (rand.event(p)) J.push_back(i);
        }
        // set empty? select random element
        if (!J.size()) J.push_back(rand.ru_int(0,ndim-1));

        // perturb decision variables
        vec1D xtest=xbest;
        for (auto k:J) {
          xtest[k]=gen_norm(xtest[k],pb[k],sigma);
          assert(xtest[k]>pb[k].xmin && xtest[k]<pb[k].xmax);
        }
        double ftest=func(xtest);

        #ifndef DDS_SIGMA_ADAPT
          if (ftest<fbest) {
            fbest=ftest;
            xbest=xtest;
          }
        #else
          if (ftest<fbest) {
            fbest=ftest;
            xbest=xtest;

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
        if (verbose) std::cout << " DDS " << std::format("{:5}",nfunc) << ": " << std::format("{:0.4f}",fbest) << " s=" << sigma << "\r";
      }
      if (verbose) std::cout << '\n';
      return {fbest,xbest};
    }
  protected:
    const int ndim;
};

#endif // DDS_H

