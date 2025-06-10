#include "cma.h"
#include "ssc.h"
#include "../common/math.h"

OptCMA::OptCMA(const CMACfg &cfg,const box_const &parambox,bool verbose)
:Opt(parambox),cfg(cfg),p(ndim),
chol(ndim),
verbose(verbose)
{
  p.sigma = cfg.sigma_init;
  p.psucc = p.p_target_succ;
}


// generate mvn-distributed variable
// xgen=x + sigma*(G*z)
auto OptCMA::generate_candidate(const vec1D &x,double sigma)
{
  vec1D xnorm=scale(x);

  vec1D z(ndim);
  for (auto &r:z)
    r = rand.r_norm();

  vec1D az=slmath::mul(chol.G,z);

  vec1D xgen(ndim);
  for (int i=0;i<ndim;i++)
  {
    #if 1
      double scale=(pb[i].xmax-pb[i].xmin)*sigma;
      double xnew=x[i] + scale*az[i];
      xgen[i]=reflect(xnew,pb[i].xmin,pb[i].xmax);
    #else
      double xnew=xnorm[i] + sigma*az[i];
      xgen[i] = reflect(xnew,0,1);
    #endif
  }
  return std::tuple{xgen,az};
}

void OptCMA::update_cov(vec2D &mcov, vec1D &pc,const vec1D &az)
{
  pc = slmath::mul_add(1.0-p.cc, pc, std::sqrt(p.cc*(2.0-p.cc)), az);
  mcov = slmath::mul_add(1.0-p.ccov, mcov, p.ccov, slmath::outer(pc,pc));
}

Opt::ppoint OptCMA::run(opt_func func,const vec1D &xstart)
{
  vec1D pc(ndim); // evolution path
  vec2D mcov(ndim,vec1D(ndim)); //covariance matrix
  for (int i=0;i<ndim;i++)
    mcov[i][i]=1.0;


  //SSC1 ssc(0.05,0.10,0.05);
  //std::cout << p.p_target_succ << ' '  << p.cp << ' ' << 1.0/p.d << '\n';
  SSC1 ssc(p.p_target_succ,p.cp,1.0/p.d);

  int nfunc=1;
  ppoint xb{func(xstart),xstart};
  if (verbose) std::cout << xb.first << '\n';

  while (nfunc < cfg.nfunc_max)
  {
    chol.Factor(mcov,0.1);

    auto [xgen,az]=generate_candidate(xb.second,p.sigma);

    double fn = func(xgen);

    double lambda=(fn<xb.first)?1.0:0.0;
    p.sigma = ssc.update(p.sigma,lambda);

    if (fn < xb.first)
    {
      xb.first = fn;
      xb.second = xgen;

      update_cov(mcov,pc,az);
    }

    nfunc++;
    if (verbose) {
      std::cout << "CMA " << std::format("{:5}",nfunc) << ": " << std::format("{:0.2f}",xb.first) << " s=" << std::format("{:0.4f}",p.sigma);
      std::cout << "\r";
    }
  }

  return xb;
}
