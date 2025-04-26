#ifndef COST_H
#define COST_H

#include "vle.h"
#include "../common/utils.h"
#include <cmath>

class CostFunction {
  public:
    CostFunction() {};
    virtual double Calc(span_ci32 buf) const =0;
    virtual ~CostFunction(){};
};

class CostL1 : public CostFunction {
  public:
    double Calc(span_ci32 buf) const override
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
    double Calc(span_ci32 buf) const override
    {
      if (buf.size()) {
        int64_t sum=0.0;
        for (const auto val:buf)
          sum+=val*val;
        return sqrt(sum/static_cast<double>(buf.size()));
      } else return 0.;
    }
};


// estimate bytes per frame with a simple golomb model
class CostGolomb : public CostFunction {
  const double alpha=0.97; // critical
  public:
    CostGolomb(){};
    double Calc(span_ci32 buf) const override
    {
      RunWeight rm(alpha);
      if (buf.size()) {
        int64_t nbits=0;
        for (const auto sval:buf) {
          const auto m=std::max(static_cast<int32_t>(rm.sum),1);
          const auto uval=MathUtils::S2U(sval);
          int q=uval/m;
          //int r=val-q*m;
          nbits+=(q+1);
          if (m>1) {
            nbits+=BitUtils::count_bits32(m);
          }
          rm.Update(uval);
        }
        return nbits/(8.);
      } else return 0;
    }
};

//#define TOTAL_SELF_INFORMATION
// entropy using order-0 markov model
class CostEntropy : public CostFunction {
  public:
    CostEntropy(){};
    double Calc(span_ci32 buf) const override
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
        std::vector<int> counts(maxval-minval+1);

        const auto cmap=[&](int32_t val) -> int& {
          return counts[val-minval];};

        for (const auto val:buf)
          ++cmap(val);

        const double invs=1.0/static_cast<double>(buf.size());
        #ifdef TOTAL_SELF_INFORMATION
          for (const auto val:buf) {
            const double p=cmap(val)*invs;
            entropy+=p*log(p);
          }
        #else
          if (counts.size() < buf.size()) { // over alphabet
            for (const auto c:counts) {
              if (c==0) continue;
              const double p=c*invs;
              entropy += c*log2(p);
            }
          } else { // over input
            for (const auto val:buf) {
              const double p=cmap(val)*invs;
              entropy+=log2(p);
            }
          }
          entropy = -entropy / 8.0;
        #endif
      }
      return entropy;
    }
};

/*class StaticBitModel {
  public:
    StaticBitModel()
    :pr(2,vec1D(PSCALE))
    {
      for (std::size_t i=1;i<PSCALE;i++)
      {
        double p1 = static_cast<double>(i)/static_cast<double>(PSCALE);
        double t0 = -std::log2(1.0-p1);
        double t1 = -std::log2(p1);
        pr[0][i]= t0;
        pr[1][i]= t1;
      }
      ResetCount();
    }
    void ResetCount(){nbits=0;};
    void EncodeBitOne(uint32_t p1,int bit)
    {
      nbits += pr[bit][p1];
    }
    auto EncodeP1_Func() {return [&](uint32_t p1,int bit) {return EncodeBitOne(p1,bit);};}; // stupid C++

  double nbits;
  vec2D pr;
};*/

class CostBitplane : public CostFunction {
 public:
  CostBitplane() {
  }
  double Calc(span_ci32 buf) const override
  {
    int numsamples=buf.size();
    std::vector<int32_t> ubuf(numsamples);
    int vmax=0;
    for (int i=0;i<numsamples;i++) {
       int val=MathUtils::S2U(buf[i]);
       if (val>vmax) vmax=val;
       ubuf[i]=val;
    }
    #if 1
    BufIO iobuf;
    RangeCoderSH rc(iobuf);
    rc.Init();
    BitplaneCoder bc_rc(MathUtils::iLog2(vmax),numsamples);
    bc_rc.Encode(rc.encode_p1,&ubuf[0]);
    rc.Stop();
    double c0=iobuf.GetBufPos();
    #else

    StaticBitModel bm;
    BitplaneCoder bc_bit(MathUtils::iLog2(vmax),numsamples);
    bc_bit.Encode(bm.EncodeP1_Func(),&ubuf[0]);

    double c0=bm.nbits/8.0;
    #endif
    return c0;
  }
};

#endif
