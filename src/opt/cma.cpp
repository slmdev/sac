#include "cma.h"

OptCMA::OptCMA(const CMACfg &cfg,const box_const &parambox,bool verbose)
:Opt(parambox),cfg(cfg),verbose(verbose)
{

}

Opt::ppoint OptCMA::run(opt_func func,const vec1D &xstart)
{
  ppoint xb;
  return xb;
}
