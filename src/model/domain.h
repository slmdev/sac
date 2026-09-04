#ifndef _DOMAIN_H
#define _DOMAIN_H

#include "model.h"

using namespace BM;

class LogDomain {
  public:
    int fmin,fmax;
    LogDomain();
    inline int Fwd(int p) const
    {
       return FwdTbl[p];
    }
    inline int Inv(int x) const
    {
       if (x<DMIN) return 1;
       else if (x>DMAX) return PSCALEm;
       else return InvTbl[x-DMIN];
    }
    void Check();
  private:
    int FwdTbl[PSCALE];
    int InvTbl[DSCALE];
};

inline LogDomain myDomain;

#endif
