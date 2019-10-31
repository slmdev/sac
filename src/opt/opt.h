#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "../global.h"
#include "../common/rand.h"
#include <functional>

// general minimization for multivariate problems using box-constraints
class Opt {
  public:
    Opt():rand(0)
    {

    }
    struct tboxconst {
      double xmin,xmax;
    };
    typedef std::pair<double,vec1D> opt_ret;
    typedef std::vector <tboxconst> param_const;
    typedef std::function<double(const vec1D &param)> tfunc;
  protected:
    // reflect xnew at boundaries
    static double Reflect(double xnew,double xmin,double xmax) {
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
    Random rand;
};

#endif // OPT_H
