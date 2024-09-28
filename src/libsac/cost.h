#ifndef COST_H
#define COST_H

#include "vle.h"
#include "../common/utils.h"
#include <cmath>

class CostFunction {
  public:
      CostFunction(){};
      virtual double Calc(const int32_t *buf,int numsamples)=0;
      virtual ~CostFunction(){};
};

class CostL1 : public CostFunction {
  public:
      double Calc(const int32_t *buf,int numsamples)
      {
        if (numsamples) {
          int64_t sum=0;
          for (int i=0;i<numsamples;i++) sum+=std::fabs(buf[i]);
          return sum/static_cast<double>(numsamples);
        } else return 0.;
      }
};

class CostRMS : public CostFunction {
  public:
      double Calc(const int32_t *buf,int numsamples)
      {
        if (numsamples) {
          int64_t sum=0.0;
          for (int i=0;i<numsamples;i++) sum+=buf[i]*buf[i];
          return sqrt(sum/static_cast<double>(numsamples));
        } else return 0.;
      }
};


// estimate number of needed bytes with a simple golomb model
// alpha paramater is critical
class CostGolomb : public CostFunction {
  const double alpha=0.97;
  const double log2=log(2.0);
  public:
      CostGolomb()
      :rm(alpha) {};
      double Calc(const int32_t *buf,int numsamples)
      {
        int64_t nbits=0;
        if (numsamples) {
          for (int i=0;i<numsamples;i++) {
            int32_t val=MathUtils::S2U(buf[i]);
            int m=std::max(static_cast<int>(rm.sum),1);

            int q=val/m;
            //int r=val-q*m;
            nbits+=(q+1);
            if (m>1) {
              int b=std::ceil(log(m)/log2);
              nbits+=b;
            };
            rm.Update(val);
          }
          return nbits/static_cast<double>(8*numsamples);
        } else return 0;
      }
  private:
    RunWeight rm;
};

// entropy using order-0 markov model
class CostEntropyO0 : public CostFunction {
  public:
    CostEntropyO0(){};
    double Calc(const int32_t *buf,int numsamples)
    {
      if (numsamples<=0) return 0.0;
      std::vector <int>counts;
      std::vector<int32_t> e(numsamples);
      int imax=0;
      for (int i=0;i<numsamples;i++) {
        e[i]=MathUtils::S2U(buf[i]);
        if (e[i]>imax) imax=e[i];
      }
      counts.resize(imax+1);
      for (int i=0;i<numsamples;i++) {
        counts[e[i]]++;
      }
      double entropy=0.0;
      for (int i=0;i<numsamples;i++) {
        double p=counts[e[i]]/double(numsamples);

        entropy+=p*log(p);
      }
      return entropy;
    }
};

// entropy using order-0 markov model
class CostEntropyO0b : public CostFunction {
  public:
    CostEntropyO0b(){};
    double Calc(const int32_t *buf,int numsamples)
    {
      if (numsamples<=0) return 0.0;
      int32_t minval = std::numeric_limits<int32_t>::max();
      int32_t maxval = std::numeric_limits<int32_t>::min();
      for (int i=0;i<numsamples;i++) {
        const int32_t val=buf[i];
        if (val>maxval) maxval=val;
        if (val<minval) minval=val;
      }
      std::vector<int> counts(maxval-minval+1);

      for (int i=0;i<numsamples;i++) {
        int32_t idx=buf[i]-minval;
        counts[idx]++;
      }
      double entropy=0.0;
      for (int i=0;i<numsamples;i++) {
        int32_t idx=buf[i]-minval;
        double p=counts[idx]/double(numsamples);

        entropy+=p*log(p);
      }
      return entropy;
    }
};


class CostBitplane : public CostFunction {
 public:
  CostBitplane() {
  }
  double Calc(const int32_t *buf,int numsamples)
  {
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
