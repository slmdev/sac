#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "../global.h"
#include "../common/rand.h"
#include <functional>

// general minimization for multivariate problems using box-constraints
class Opt {
  struct tboxconst {
    double xmin,xmax;
  };
  public:
    inline static const std::string SLOPT_VERSION="0.1.1";
    typedef std::pair<double,vec1D> opt_ret;
    typedef std::vector <tboxconst> box_const;
    typedef std::function<double(const vec1D &param)> opt_func;
    Opt()
    :rand(0)
    {

    };
    Opt(const box_const &parambox)
    :rand(0),pb(parambox)
    {

    };
  protected:
    vec1D scale(const vec1D &x) {
      vec1D v_out(x.size());
      for (size_t i=0;i<x.size();i++)
        v_out[i] = (x[i]-pb[i].xmin) / (pb[i].xmax-pb[i].xmin);
      return v_out;
    }
    vec1D unscale(const vec1D &x)
    {
      vec1D v_out(x.size());
      for (size_t i=0;i<x.size();i++)
        v_out[i] = (x[i] * (pb[i].xmax-pb[i].xmin)) + pb[i].xmin;
      return v_out;
    }
    // reflect xnew at boundaries
    double reflect(double xnew,double xmin,double xmax) {
      if (xnew<xmin) {
        xnew=xmin+(xmin-xnew);
        if (xnew>xmax) xnew=xmin;
      }
      if (xnew>xmax) {
        xnew=xmax-(xnew-xmax);
        if (xnew<xmin) xnew=xmax;
      }
      return xnew;
    }
    double reset(double xnew,double xmin,double xmax) {
      if (xnew < xmin)
        return xmin;
      else if (xnew > xmax)
        return xmax;
      else
        return xnew;
    }
    Random rand;
    const box_const pb;
};

#endif
