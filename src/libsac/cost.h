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

class CostMeanRMS : public CostFunction {
  public:
      double Calc(const int32_t *buf,int numsamples)
      {
        if (numsamples) {
          double pow=0.0;
          for (int i=0;i<numsamples;i++) pow+=std::fabs(buf[i]);
          return pow/static_cast<double>(numsamples);
        } else return 0.;
      }
};

// estimate number of needed bits with a simple golomb model
class CostGolomb : public CostFunction {
  public:
      CostGolomb():mean_err(0.98){};
      double Calc(const int32_t *buf,int numsamples)
      {
        const double log2=log(2.0);
        double nbits=0;
        if (numsamples) {
          for (int i=0;i<numsamples;i++) {
            int32_t val=MathUtils::S2U(buf[i]);
            int m=(std::max)(static_cast<int>(mean_err.Get()),1);
            int q=val/m;
            //int r=val-q*m;
            nbits+=(q+1);
            if (m>1) {
              int b=ceil(log(m)/log2);
              nbits+=b;
            }
            mean_err.Update(val);
          }
          return nbits/(double)numsamples;
        } else return 0;
      }
  private:
    RunExp mean_err;
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

#endif // COST_H
