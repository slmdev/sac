#ifndef DDS_H
#define DDS_H

#include "opt.h"
#include <cassert>

// Dynamical dimensioned search algorithm for computationally efficient watershed model calibration
// Tolson, Shoemaker 2007
class DDS : public Opt {
  public:
    DDS(int nparam)
    :n(nparam)
    {
      defparam.resize(n);
      for (int i=0;i<n;i++) {defparam[i].xmin=0;defparam[i].xmax=1;};
    }
    DDS(const param_const &defparam)
    :n(defparam.size()),defparam(defparam)
    {

    }
    opt_ret Run(double r,int nfunc_max,const vec1D &param,tfunc func)
    {
      int nfunc=1;
      double fbest=func(param); // eval at initial solution
      //std::cout << fbest << '\n';
      vec1D xbest=param;

      while (nfunc<nfunc_max) {
        std::vector <int>J; // select J of D variables
        double p=1.0-log(nfunc)/log(nfunc_max);
        for (int i=0;i<n;i++) {
          if (rand.event(p)) J.push_back(i);
        }
        // set empty? select random element
        if (!J.size()) J.push_back(rand.ru_int(0,n-1));

        // perturb decision variables
        vec1D xtest=xbest;
        for (auto k:J) {
          xtest[k]=GenerateCandidate(xtest[k],defparam[k].xmin,defparam[k].xmax,r);
          assert(xtest[k]>defparam[k].xmin && xtest[k]<defparam[k].xmax);
        }
        double ftest=func(xtest);
        if (ftest<fbest) {
          fbest=ftest;
          xbest=xtest;
          /*std::cout << fbest << '\n';
          for (auto &x:xbest) std::cout << x << ' ';
          std::cout << '\n';*/
        }
        nfunc++;
        std::cout << "dds: " << nfunc << '\r';
      }
      return opt_ret{fbest,xbest};
    }
  protected:
    double GenerateCandidate(double x,double xmin,double xmax,double r)
    {
      double sigma=r*(xmax-xmin);
      double xnew=x+sigma*rand.r_norm(0,1);
      return Reflect(xnew,xmin,xmax);
    }
    int n;
    param_const defparam;
};

#endif // DDS_H
