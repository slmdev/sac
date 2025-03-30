#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "../global.h"
#include "../common/rand.h"
#include <functional>
#include <cassert>

// general minimization for multivariate problems using box-constraints
class Opt {
  public:
    struct tboxconst {
      double xmin,xmax;
    };
    inline static const std::string SLOPT_VERSION="0.3.0";
    typedef std::pair<double,vec1D> ppoint;
    typedef std::vector<ppoint> opt_points;
    typedef std::vector <tboxconst> box_const;
    typedef std::function<double(const vec1D &param)> opt_func;

    Opt(const box_const &parambox);

  protected:
    std::size_t eval_points_mt(opt_func func,std::span<ppoint> ps);

    // scale to [0,1]
    vec1D scale(const vec1D &x);
    vec1D unscale(const vec1D &x);
    double unscale(double r,const tboxconst &box);

    // generate random normal distributed sample around x with sigma r
    double gen_norm(const double x,const tboxconst &box,const double r);

    vec1D gen_norm_samples(const vec1D &xb,double r);
    vec1D gen_uniform_samples(const vec1D &xb, double r);

    vec1D gen_uniform_sample();
    // reflect xnew at boundaries
    double reflect(double xnew,double xmin,double xmax);
    double reset(double xnew,double xmin,double xmax);

    Random rand;
    const box_const pb;
    const int ndim;
};

#endif
