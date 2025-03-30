#include "opt.h"
#include <future>

Opt::Opt(const box_const &parambox)
:rand(0),pb(parambox),ndim(parambox.size())
{

};

// evaluate span of points with ps.size() parallel threads
std::size_t Opt::eval_points_mt(opt_func func,std::span<ppoint> ps)
{
  std::vector <std::future<double>> threads;
  threads.reserve(ps.size());

  for (std::size_t i=0;i<ps.size();i++) {
    auto arg=ps[i].second;
    threads.emplace_back(std::async(std::launch::async, [func, arg]() {
      return func(arg);
    }));
  }

  if (threads.size()!=ps.size())
    std::cerr << "  warning: eval_points_mt: thread count too low\n";

  std::vector<double>r1(ps.size());
  for (std::size_t i=0;i<threads.size();i++) {
    r1[i]=threads[i].get();
    if (std::isnan(r1[i])) std::cerr << " warning: nan in eval_points_mt\n";
  }

  // transfer savely to ps
  for (std::size_t i=0;i<ps.size();i++)
    ps[i].first=r1[i];

  #if 0 // check for thread safety
    std::vector<double> r2(ps.size());
    for (std::size_t i=0;i<ps.size();i++)
      r2[i] = func(ps[i].second);

    for (std::size_t i=0;i<ps.size();i++)
      if (r1[i] != r2[i])
        std::cerr << "  warning: mt res (" << i << "): " << r1[i] << ' ' << r2[i] << '\n';
  #endif

  return threads.size();
}

vec1D Opt::scale(const vec1D &x) {
  vec1D v_out(x.size());
  for (size_t i=0;i<x.size();i++)
    v_out[i] = (x[i]-pb[i].xmin) / (pb[i].xmax-pb[i].xmin);
  return v_out;
}
vec1D Opt::unscale(const vec1D &x)
{
  vec1D v_out(x.size());
  for (size_t i=0;i<x.size();i++)
    v_out[i] = (x[i] * (pb[i].xmax-pb[i].xmin)) + pb[i].xmin;
  return v_out;
}

// generate random normal distributed sample around x with sigma r
double Opt::gen_norm(const double x,const tboxconst &box,const double r)
{
  double sigma=r*(box.xmax-box.xmin);
  double xnew=x+sigma*rand.r_norm(0,1);
  return reflect(xnew,box.xmin,box.xmax);
}
double Opt::unscale(double r,const tboxconst &box)
{
  return (r*(box.xmax-box.xmin)) + box.xmin;
}

vec1D Opt::gen_norm_samples(const vec1D &xb,double radius)
{
  assert(pb.size()==xb.size());
  const int n=pb.size();
  vec1D v_out(n);
  for (int i=0;i<n;i++)
    v_out[i] = gen_norm(xb[i],pb[i],radius);
  return v_out;
}
vec1D Opt::gen_uniform_samples(const vec1D &xb, double radius)
{
  assert(pb.size()==xb.size());
  const int n=pb.size();
  vec1D v_out(n);
  for (int i=0;i<n;i++) {
    double s=radius*(pb[i].xmax-pb[i].xmin);
    double rnd=rand.r_int(-1,+1);
    double xn=xb[i] + rnd*s;
    v_out[i] = reflect(xn,pb[i].xmin,pb[i].xmax);
  }
  return v_out;
}

vec1D Opt::gen_uniform_sample()
{
  int n=pb.size();
  vec1D v_out(n);
  for (int i=0;i<n;i++) {
    v_out[i] = unscale(rand.r_01closed(),pb[i]); //sample from [0,1] and scale
  }
  return v_out;
}

// reflect xnew at boundaries
double Opt::reflect(double xnew,double xmin,double xmax) {
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

double Opt::reset(double xnew,double xmin,double xmax) {
  if (xnew < xmin)
    return xmin;
  else if (xnew > xmax)
    return xmax;
  else
    return xnew;
}

