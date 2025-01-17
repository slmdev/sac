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
    inline static const std::string SLOPT_VERSION="0.2.0";
    typedef std::pair<double,vec1D> ppoint;
    typedef std::vector<ppoint> opt_points;
    typedef std::vector <tboxconst> box_const;
    typedef std::function<double(const vec1D &param)> opt_func;

    Opt(const box_const &parambox)
    :rand(0),pb(parambox)
    {

    };
  protected:
    // scale to [0,1]
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
    // generate random normal distributed sample around x with sigma r
    double gen_norm(const double x,const tboxconst &box,const double r)
    {
      double sigma=r*(box.xmax-box.xmin);
      double xnew=x+sigma*rand.r_norm(0,1);
      return reflect(xnew,box.xmin,box.xmax);
    }
    double unscale(double r,const tboxconst &box)
    {
      return (r*(box.xmax-box.xmin)) + box.xmin;
    }

    vec1D gen_norm_sample(const vec1D &xb,double r)
    {
      assert(pb.size()==xb.size());
      const int n=pb.size();
      vec1D v_out(n);
      for (int i=0;i<n;i++)
        v_out[i] = gen_norm(xb[i],pb[i],r);
      return v_out;
    }
    vec1D gen_uniform_sample(const vec1D &xb, double r)
    {
      assert(pb.size()==xb.size());
      const int n=pb.size();
      vec1D v_out(n);
      for (int i=0;i<n;i++) {
        double s=r*(pb[i].xmax-pb[i].xmin);
        double r=rand.r_int(-1,+1);
        double xn=xb[i] + r*s;
        v_out[i] = reflect(xn,pb[i].xmin,pb[i].xmax);
      }
      return v_out;
    }

    vec1D gen_uniform_sample()
    {
      int n=pb.size();
      vec1D v_out(n);
      for (int i=0;i<n;i++) {
        v_out[i] = unscale(rand.r_01closed(),pb[i]); //sample from [0,1] and scale
      }
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
