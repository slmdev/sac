#ifndef MODEL_H
#define MODEL_H

#include <algorithm>

// probability precision

namespace BM {
  //probability
  inline constexpr int PBITS  = 15;
  inline constexpr int PSCALE = 1<<PBITS;
  inline constexpr int PSCALEh= PSCALE>>1;
  inline constexpr int PSCALEm= PSCALE-1;

  // mixer weights
  inline constexpr int WBITS  = 16; //fractional bits
  inline constexpr int WSCALE = 1<<WBITS;
  inline constexpr int WRANGE = 1<<(WBITS+3);

  //mixer learning rate
  inline constexpr int RBITS  = 16;
  inline constexpr int RSCALE = 1 << RBITS;

  inline constexpr int MixerRate(double eta) {return int(eta*RSCALE+0.5);};

  //logit domain
  inline constexpr int LBITS = 8; //fractional bits
  inline constexpr int LSCALE = 1<<LBITS;

  inline constexpr int LOGITMAX = 8; //saturation point
  inline constexpr int DMAX = LOGITMAX*LSCALE-1;
  inline constexpr int DMIN = -DMAX;
  inline constexpr int DSCALE = DMAX-DMIN+1;

  //total mixer shift
  inline constexpr int USHIFT = LBITS + PBITS - WBITS + RBITS;
}

#endif
