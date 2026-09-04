#include <cstdio>
#include <cmath>

#include "domain.h"

LogDomain::LogDomain()
{
  for (int i=0;i<PSCALE;i++)
  {
    double f=std::max(i,1)/(double)PSCALE;
    double q=std::log(f / (1.0-f))*LSCALE;
    FwdTbl[i]=(int)std::round(q);
  };
  fmin=FwdTbl[0];
  fmax=FwdTbl[PSCALE-1];
  for (int i=DMIN;i<=DMAX;i++)
  {
    double q=PSCALE/(1.0+exp(-double(i)/double(LSCALE)));
    InvTbl[i-DMIN]=(int)std::round(q);
  };
}

void LogDomain::Check()
{
  int sum=0;
  for (int i=0;i<PSCALE;i++)
  {
    int p=Inv(Fwd(i));
    sum+=(p-i)*(p-i);
  }
  printf(" mse: %0.4f\n",double(sum)/double(PSCALE));
}
