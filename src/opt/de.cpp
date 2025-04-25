#include <cassert>
#include "de.h"
#include "../common/utils.h"

OptDE::OptDE(const DECfg &cfg,const box_const &parambox,bool verbose)
:Opt(parambox),cfg(cfg),verbose(verbose)
{
}

// select k unique elements from 0..n-1 except ie
std::vector<int> OptDE::select_k_unique_except(int n,int ie,int k)
{
  std::vector<int>r;
  if (k>=n-1) return r;

  std::vector<int> e(n);
  std::iota(std::begin(e),std::end(e),0);
  std::erase(e,ie);

  for (int i=0;i<k;i++)
  {
    int idx = rand.ru_int(0,e.size()-1);
    int val = e[idx];
    r.push_back(val);
    std::erase(e,val);
  }
  return r;
}

auto OptDE::generate_candidate(const opt_points &pop,const vec1D &xbest,int iagent,double mCR,double mF)
{
  const double tCR = gen_CR(mCR);
  const double tF  = gen_F(mF);
  const int R = rand.ru_int(0,ndim-1);

  auto gp = [&](int i) -> auto& {return pop[i].second;};

  // select k distinct elements from pop depending on mut_method
  auto v = select_k_unique_except(pop.size(),iagent,MutVals[cfg.mut_method]);

  // mutation
  vec1D xm;
  if (cfg.mut_method==BEST1BIN) {
    xm = mut_1bin(xbest,gp(v[0]),gp(v[1]),mF);
  } else if (cfg.mut_method==RAND1BIN) {
    xm = mut_1bin(gp(v[0]),gp(v[1]),gp(v[2]),mF);
  } else if (cfg.mut_method==CUR1BEST) {
    xm = mut_curbest(xbest,gp(iagent),gp(v[0]),gp(v[1]),mF);
  } else if (cfg.mut_method==CURPBEST) {
    int np = std::min(cfg.npbest,static_cast<int>(pop.size())-1); // pop.size() can by smaller than cfg.NP
    int xp = rand.ru_int(0,np);
    xm = mut_curbest(gp(xp),gp(iagent),gp(v[0]),gp(v[1]),mF);
  }

  // cross-over
  vec1D xtrial(ndim);
  const ppoint &xi=pop[iagent];
  for (int i=0;i<ndim;i++)
  {
    if (rand.event(mCR) || (i==R)) xtrial[i] = xm[i];
    else xtrial[i] = xi.second[i];
  }

  return std::tuple{xtrial,tCR,tF};
}

void OptDE::print_status(std::size_t nfunc,double fx,double mCR,double mF)
{
  std::cout << "DE " << std::format("{:5}",nfunc) << ": " << std::format("{:0.2f}",fx);
  std::cout << " (mCR=" << std::format("{:0.4f}",mCR) << " mF=" << std::format("{:0.4f}",mF);
  std::cout << ")\r";
}

OptDE::ppoint OptDE::run(opt_func func,const vec1D &xstart)
{
  assert(pb.size()==xstart.size());

  std::size_t nfunc=1;
  ppoint xb{func(xstart),xstart}; // eval at initial solution
  if (verbose) std::cout << xb.first << '\n';

  opt_points pop(cfg.NP); // population
  pop[0] = xb;

  std::span<ppoint> pop_span(pop.begin()+1,pop.end());
  // generate random population
  for (auto &x : pop_span) {
    vec1D xt;
    if (cfg.init_method == INIT_UNIV)
      xt=gen_uniform_samples(xb.second,cfg.sigma_init);
    else if (cfg.init_method == INIT_NORM)
      xt=gen_norm_samples(xb.second,cfg.sigma_init);

    x.second = xt;
  }

  // eval inital population in parallel (excluding first)
  nfunc +=eval_pop_pool(func,pop_span,cfg.num_threads);
  // update best sample
  for (const auto &x:pop_span)
    if (x.first < xb.first)
      xb = x;

  double mCR = cfg.CR;
  double mF = cfg.F;

  if (verbose) std::cout << "DE init pop " << pop.size() << " (mt=" << cfg.num_threads << "): s=" << cfg.sigma_init << ": " << xb.first << "\n";

  // trial agents
  opt_points gen_pop;
  std::vector<std::pair<double,double>> gen_mut;

  while (nfunc<cfg.nfunc_max)
  {
    if (cfg.mut_method==CURPBEST) { // sort by function value
      std::sort(begin(pop),end(pop),
        [](const auto &a,const auto &b){return a.first < b.first;});
    }

    int num_agents = std::min(cfg.nfunc_max-nfunc,pop.size());

    // trial agents
    gen_mut.resize(num_agents);
    gen_pop.resize(num_agents);

    for (int iagent=0;iagent<num_agents;iagent++)
    {
      auto [xtrial, tCR, tF] = generate_candidate(pop,xb.second,iagent,mCR,mF);
      gen_mut[iagent]={tCR,tF}; // save (random) mutation params
      gen_pop[iagent].second = xtrial;
    }

    // evaluate trial population
    nfunc+=eval_pop_pool(func,gen_pop,cfg.num_threads);

    // greedy selection
    std::vector<double>CR_succ;
    std::vector<double>F_succ;
    for (int iagent=0;iagent<num_agents;iagent++)
      if (gen_pop[iagent].first < pop[iagent].first)
      {
        pop[iagent] = gen_pop[iagent]; // replace
        CR_succ.push_back(gen_mut[iagent].first);
        F_succ.push_back(gen_mut[iagent].second);

        // replace best vector
        if (pop[iagent].first < xb.first)
          xb = pop[iagent];
      }

    if (verbose) print_status(nfunc,xb.first,mCR,mF);

    if (nfunc >= cfg.nfunc_max) break;

    mCR = (1.0-cfg.c)*mCR + cfg.c*MathUtils::mean(CR_succ);
    mF  = (1.0-cfg.c)*mF  + cfg.c*MathUtils::mean(F_succ);
  }
  if (verbose) std::cout << '\n';
  return xb;
}

vec1D OptDE::mut_1bin(const vec1D &xb,const vec1D &x1,const vec1D &x2,double F)
{
  vec1D xm(ndim);
  for (int i=0;i<ndim;i++) {
    double y=xb[i] + F*(x1[i] - x2[i]);
    xm[i] = reflect(y, pb[i].xmin,pb[i].xmax);
  }
  return xm;
}

vec1D OptDE::mut_curbest(const vec1D &xbest,const vec1D &xb,const vec1D &x1,const vec1D &x2,double F)
{
  vec1D xm(ndim);
  for (int i=0;i<ndim;i++) {
    double y=xb[i] + F*(xbest[i] - xb[i]) + F*(x1[i] - x2[i]);
    xm[i] = reflect(y, pb[i].xmin,pb[i].xmax);
  }
  return xm;
}
