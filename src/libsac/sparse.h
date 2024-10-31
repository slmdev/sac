#ifndef SPARSEPCM_H
#define SPARSEPCM_H

class SimplePred {
  public:
    SimplePred()
    :lb(0)
    {

    }
    double Predict()
    {
      return lb;
    }
    void Update(int32_t val)
    {
      lb = val;
    }
  protected:
    int32_t lb;
};

class SparsePCM {
  const double cost_pow=1;
  public:
    SparsePCM()
    :minval(0),maxval(0),fraction_used(0.),fraction_cost(0.)
    {
    };
    void Analyse(std::span<const int32_t> buf)
    {
      minval = std::numeric_limits<int32_t>::max();
      maxval = std::numeric_limits<int32_t>::min();
      for (auto val : buf) {
        if (val>maxval) maxval=val;
        if (val<minval) minval=val;
      }
      used.resize(maxval-minval+1);

      for (auto val : buf) used[val-minval] = 1;
      int sum=0;
      for (auto u : used) sum+=u;
      fraction_used = used.size()>0?(sum/static_cast<double>(used.size()))*100.:0.0;


      // calc cost
      //SimplePred pred;
      double sum0=0,sum1=0;
      for (auto val : buf) {
        //int32_t p=std::clamp((int)std::round(pred.Predict()),minval,maxval);
        int32_t e0=val;
        int32_t e1=map_val(e0);

        sum0+=pow(std::fabs(e0),cost_pow);
        sum1+=pow(std::fabs(e1),cost_pow);

        //pred.Update(val);
      }
      fraction_cost=sum1>0?sum0/static_cast<double>(sum1):0;
    }
    int map_val(const int32_t val,const int32_t p=0)
    {
      if (val==0) return 0;
      const int sgn=MathUtils::sgn(val);

      const int pidx=p-minval;
      int mres=0;
      for (int i=1;i<=sgn*val;i++)
      {
         const int idx=pidx+(sgn*i);
         /*if (idx<0 || (unsigned)idx>=used.size())
         {
           std::cout << p << ' ' << val << ' ' << p+(sgn*i) << ' ' << minval << ' ' << maxval << '\n';
         } else */
         if (used[idx]) ++mres;
      }
      return sgn*mres;
    }
    int32_t minval,maxval;
    double fraction_used,fraction_cost;
  protected:
    std::vector<int> used;
};


#endif // SPARSEPCM_H
