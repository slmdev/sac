#ifndef COST_H
#define COST_H

#include "vle.h"
#include "../common/utils.h"
#include <cmath>

typedef std::span<int32_t> span_i32;

class CostFunction {
  public:
      CostFunction() {};
      virtual double Calc(const span_i32 &buf) const =0;
      virtual ~CostFunction(){};
};

class CostL1 : public CostFunction {
  public:
      double Calc(const span_i32 &buf) const
      {
        if (buf.size()) {
          int64_t sum=0;
          for (const auto val:buf)
            sum+=std::fabs(val);
          return sum/static_cast<double>(buf.size());
        } else return 0.;
      }
};

class CostRMS : public CostFunction {
  public:
      double Calc(const span_i32 &buf) const
      {
        if (buf.size()) {
          int64_t sum=0.0;
          for (const auto val:buf)
            sum+=val*val;
          return sqrt(sum/static_cast<double>(buf.size()));
        } else return 0.;
      }
};


// estimate number of needed bytes with a simple golomb model
// alpha paramater is critical
class CostGolomb : public CostFunction {
  const double alpha=0.97;
  const double log2=log(2.0);
  public:
      CostGolomb(){};
      double Calc(const span_i32 &buf) const
      {
        RunWeight rm(alpha);
        int64_t nbits=0;
        if (buf.size()) {
          for (const auto sval:buf) {
            const auto m=std::max(static_cast<int32_t>(rm.sum),1);
            const auto uval=MathUtils::S2U(sval);
            int q=uval/m;
            //int r=val-q*m;
            nbits+=(q+1);
            if (m>1) nbits+=std::ceil(log(m)/log2);
            rm.Update(uval);
          }
          return nbits/static_cast<double>(8*buf.size());
        } else return 0;
      }
};

// entropy using order-0 markov model
class CostEntropy : public CostFunction {
  public:
    CostEntropy(){};
    double Calc(const span_i32 &buf) const
    {
      double entropy=0.0;
      if (buf.size())
      {
        int32_t minval = std::numeric_limits<int32_t>::max();
        int32_t maxval = std::numeric_limits<int32_t>::min();
        for (const auto val:buf) {
          if (val>maxval) maxval=val;
          if (val<minval) minval=val;
        }
        auto vmap=[&](int32_t val) {return val-minval;};

        std::vector<int> counts(maxval-minval+1);
        for (const auto val:buf)
          counts[vmap(val)]++;

        const double invs=1.0/static_cast<double>(buf.size());
        for (const auto val:buf) {
          const double p=counts[vmap(val)]*invs;
          entropy+=p*log(p);
        }
      }
      return entropy;
    }
};


class CostBitplane : public CostFunction {
 public:
  CostBitplane() {
  }
  double Calc(const span_i32 &buf) const
  {
    int numsamples=buf.size();
    std::vector<int32_t> ubuf(numsamples);
    int vmax=0;
    for (int i=0;i<numsamples;i++) {
       int val=MathUtils::S2U(buf[i]);
       if (val>vmax) vmax=val;
       ubuf[i]=val;
    }
    BufIO iobuf;
    RangeCoderSH rc(iobuf);
    rc.Init();
    BitplaneCoder bc(rc,MathUtils::iLog2(vmax),numsamples);
    bc.Encode(&ubuf[0]);
    rc.Stop();
    return iobuf.GetBufPos();
  }
};

#endif
