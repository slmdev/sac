#include "opt.h"
#include <cassert>

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

  for (std::size_t i=0;i<threads.size();i++) {
    ps[i].first=threads[i].get();
    if (std::isnan(ps[i].first))
      std::cerr << " warning: nan in eval_points_mt\n";
  }

  #if 0 // check for thread safety
    std::vector<double> rt(ps.size());
    for (std::size_t i=0;i<ps.size();i++)
      rt[i] = func(ps[i].second);

    for (std::size_t i=0;i<ps.size();i++)
      if (ps[i].first != rt[i])
        std::cerr << "  warning: mt res (" << i << "): " << ps[i].first << ' ' << rt[i] << '\n';
  #endif

  return threads.size();
}

// evaluate population parallel in rounds of num_threads
// all threads should have the same work load
std::size_t Opt::eval_pop(opt_func func,std::span<ppoint> pop,std::size_t num_threads)
{
  std::size_t n=0;
  while (n < pop.size())
  {
    std::size_t start=n;
    std::size_t ende=std::min(pop.size(),n+num_threads);
    std::size_t k=eval_points_mt(func,std::span{begin(pop)+start,begin(pop) + ende});

    n+=k;
  }
  return n;
}

// evaluate population with a simple thread-pool using a shared atomic counter
// more efficient if work load is different per thread
std::size_t Opt::eval_pop_pool(opt_func func,std::span<ppoint> pop,std::size_t num_threads)
{
  // shared atomic counter among threads
  std::atomic<std::size_t> index{0};

  auto worker = [&]() {
    while (true) {
      // return counter, inc after - only need atomic uniqueness not synchronization visibility
      std::size_t i = index.fetch_add(1,std::memory_order_relaxed);
      if (i >= pop.size()) break; // index oob - nothing to do

      double result = func(pop[i].second);
      pop[i].first = result;

      if (std::isnan(result)) {
          std::cerr << " warning: nan in eval_pop\n";
      }
    }
  };

  // launch workers
  std::vector<std::thread> threads;
  for (std::size_t i=0;i<std::min(pop.size(),num_threads);i++) {
    threads.emplace_back(worker);
  }

  for (auto &t:threads) {
    t.join();
  }

  return pop.size();
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
  double xnew=x+sigma*rand.r_norm();
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

