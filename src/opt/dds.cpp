#include <format>
#include "dds.h"
#include "ssc.h"


OptDDS::OptDDS(const DDSCfg &cfg,const box_const &parambox,bool verbose)
:Opt(parambox),cfg(cfg),
verbose(verbose)
{
}

vec1D OptDDS::generate_candidate(const vec1D &x,int nfunc,double sigma)
{
  std::vector <int>J; // select J of D variables
  double p=1.0-log(nfunc)/log(cfg.nfunc_max);

  for (int i=0;i<ndim;i++) {
    if (rand.event(p)) J.push_back(i);
  }
  // set empty? select random element
  if (!J.size()) J.push_back(rand.ru_int(0,ndim-1));

  // perturb decision variables
  vec1D xtest=x;
  for (auto k:J) {
    xtest[k]=gen_norm(x[k],pb[k],sigma);
    assert(xtest[k]>=pb[k].xmin && xtest[k]<=pb[k].xmax);
  }
  return xtest;
}

// sequential single threaded
Opt::ppoint OptDDS::run_single(opt_func func,const vec1D &xstart)
{
  int nfunc=1;
  ppoint xb{func(xstart),xstart};
  if (verbose) std::cout << xb.first << '\n';

  double sigma=cfg.sigma_init;

  // step size control
  SSC0 ssc(cfg.c_succ_max,cfg.c_fail_max);

  while (nfunc<cfg.nfunc_max) {
    ppoint x_gen;
    x_gen.second=generate_candidate(xb.second,nfunc,sigma);
    x_gen.first=func(x_gen.second);
    nfunc++;

    double lambda=0.0;
    if (x_gen.first<xb.first) {
      xb = x_gen;
      lambda=1.0;
    }
    sigma = ssc.update(sigma,lambda);

    if (verbose) std::cout << " DDS " << std::format("{:5}",nfunc) << ": " << std::format("{:0.4f}",xb.first) << " s=" << sigma << "\r";
  }
  return xb;
}

// multi-threaded variant
Opt::ppoint OptDDS::run_mt(opt_func func,const vec1D &xstart)
{

  ppoint xb{func(xstart),xstart}; // eval at initial solution

  if (verbose) std::cout << xb.first << '\n';

  double sigma=cfg.sigma_init;

  // step size control
  SSC1 ssc(0.05,0.10,0.05);

  int nfunc=1;
  while (nfunc<cfg.nfunc_max) {
    const int nthreads = std::min(cfg.nfunc_max-nfunc,cfg.num_threads);

    // generate candidates around current xbest
    opt_points x_gen(nthreads);
    for (int i=0;i<nthreads;i++) {
      x_gen[i].second=generate_candidate(xb.second,nfunc,sigma);
      nfunc++;
    };

    eval_points_mt(func,x_gen);

    // select
    ppoint xb_old=xb;
    int nsucc=0;
    for (const auto &xg : x_gen)
      if (xg.first<xb_old.first)  {
        nsucc++; // count as success, if better than parent
        if (xg.first < xb.first) // replace overall best
          xb = xg;
      }
    double lambda=nsucc/static_cast<double>(nthreads);
    sigma=ssc.update(sigma,lambda);

    if (verbose) {
        std::cout << " DDS mt=" << nthreads << ": " << std::format("{:5}",nfunc) << ": " << std::format("{:0.2f}",xb.first);
        std::cout << " s=" << std::format("{:0.3f}",sigma) << ", p_succ=" << std::format("{:0.4f}",ssc.p_succ) << "\r";
    }
  }
  return xb;
}

OptDDS::ppoint OptDDS::run(opt_func func,const vec1D &xstart)
{
  assert(pb.size()==xstart.size());

  ppoint pbest;
  if (cfg.num_threads<=0) pbest=run_single(func,xstart);
  else pbest=run_mt(func,xstart);

  if (verbose) std::cout << '\n';
  return pbest;
}
