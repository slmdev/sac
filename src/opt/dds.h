#ifndef DDS_H
#define DDS_H

#include "opt.h"
#include <cassert>

// Dynamical dimensioned search algorithm for computationally efficient watershed model calibration
// Tolson, Shoemaker 2007
class DDS : public Opt {
  public:

    DDS(const box_const &parambox)
    :Opt(parambox)
    {
    }
    double cradius(double rmin, double rmax, int n, int nmax)
    {
        double dx = nmax;
        double dy = rmin-rmax;
        return n*(dy/dx)+rmax;
    }
    opt_ret run(opt_func func,const vec1D &xstart,int nfunc_max,double radius=0.1)
    {
      assert(pb.size()==xstart.size());

      int nfunc=1;
      int n=xstart.size();
      vec1D xbest=xstart;
      double fbest=func(xstart); // eval at initial solution

      //const double phi=std::min(20.0/double(n),1.0);

      while (nfunc<nfunc_max) {
        std::vector <int>J; // select J of D variables
        double p=(1.0-log(nfunc)/log(nfunc_max));

        for (int i=0;i<n;i++) {
          if (rand.event(p)) J.push_back(i);
        }
        // set empty? select random element
        if (!J.size()) J.push_back(rand.ru_int(0,n-1));

        // perturb decision variables
        vec1D xtest=xbest;
        for (auto k:J) {
          xtest[k]=GenerateCandidate(xtest[k],pb[k].xmin,pb[k].xmax,cradius(0.05,0.20,nfunc,nfunc_max));
          assert(xtest[k]>pb[k].xmin && xtest[k]<pb[k].xmax);
        }
        double ftest=func(xtest);
        if (ftest<fbest) {
          fbest=ftest;
          xbest=xtest;
        }
        nfunc++;
      }
      return {fbest,xbest};
    }
  protected:
    double GenerateCandidate(double x,double xmin,double xmax,double r)
    {
      double sigma=r*(xmax-xmin);
      double xnew=x+sigma*rand.r_norm(0,1);
      return reflect(xnew,xmin,xmax);
    }
};

#endif // DDS_H

