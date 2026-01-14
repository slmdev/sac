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

// rank-mapping errors
// counting how many used symbols exist between the prediction and the target: error+prediction
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

      used.assign(maxval-minval+1,0);

      for (auto val : buf) used[val-minval] = 1;
      int sum=std::accumulate(begin(used),end(used),0);
      fraction_used = used.size()>0?(sum/static_cast<double>(used.size()))*100.:0.0;

      //create prefix sum, count used values before "i"
      prefix.assign(used.size()+1,0);
      for (size_t i=0;i<used.size();++i) {
        prefix[i+1]=prefix[i]+used[i];
      }

      // calc cost
      //SimplePred pred;
      double sum0=0,sum1=0;
      for (auto val : buf) {
        //int32_t p=std::clamp((int)std::round(pred.Predict()),minval,maxval);
        int32_t e0=val;
        //int32_t e1=val2rank(e0);
        int32_t e1=val2rank_fast(e0);
        //if (e1!=e2)
        //  std::cerr << "  warning: " << e1 << ' ' << e2 << '\n';

        sum0+=pow(std::fabs(e0),cost_pow);
        sum1+=pow(std::fabs(e1),cost_pow);

        //pred.Update(val);
      }
      fraction_cost=sum1>0?sum0/static_cast<double>(sum1):0;
    }

    // map error "val" (relative to prediction p) to the rank distance among only used symbols
    // mapping is bijective over valid target values val+p in [minval,maxval]
    // multiple p outside the support are semantically identical with respect to ranking inside support
    int val2rank_fast(const int32_t val,const int32_t p=0)
    {
        const int tidx = (val+p) - minval; // target index
        if (tidx < 0 || tidx>=(int)used.size() || used[tidx]==0)
          throw std::runtime_error("map_val: target index outside support");

        if (val==0) return 0;

        const int N=(int)used.size();
        const int pidx = p - minval; // pivot/prediction index

        if (val>0) {
          // range (p,target) -> prefix[tidx+1] - prefix[clamp(pidx+1)]
          int pidx_clamp = std::clamp(pidx+1,0,N); // no violation of reversibility
          return prefix[tidx+1] - prefix[pidx_clamp];
        } else {
          // range [target,p) -> prefix[clamp(pidx)] - prefix[tidx]
          int pidx_clamp = std::clamp(pidx,0,N);
          //return negative rank
          return prefix[tidx] - prefix[pidx_clamp];
        }
    }
    //possible oob
    int val2rank(const int32_t val,const int32_t p=0)
    {
      if (val==0) return 0;
      const int sgn=MathUtils::sgn(val);

      const int pidx=p-minval;
      int mres=0;
      if (val>0) {
        for (int i=pidx+1;i<=pidx+val;i++)
          mres+=used[i];
          //if (used[i]) ++mres;
      } else {
        for (int i=pidx-1;i>=pidx+val;i--)
          mres+=used[i];
         //if (used[i])  ++mres;
      }
      return sgn*mres;
    }
    int32_t minval,maxval;
    double fraction_used,fraction_cost;
  protected:
    std::vector<int> used,prefix;
};


#endif // SPARSEPCM_H
