#ifndef UTILS_H
#define UTILS_H

#include "../global.h"
#include "math.h"

#include <algorithm>
#include <numeric>
#include <cstring>
#include <cmath>
#include <immintrin.h>

// running exponential smoothing
// sum=alpha*sum+(1.0-alpha)*val, where 1/(1-alpha) is the mean number of samples considered
class RunExp {
  public:
      RunExp(double alpha):sum(0.0),alpha(alpha){};
      RunExp(double alpha,double sum):sum(sum),alpha(alpha){};
      inline void Update(double val) {
        sum=alpha*sum+(1.-alpha)*val;
      }
    double sum;
  private:
    double alpha;
};

// running weighted sum: sum_{i=0}^n alpha^(n-i) val
class RunWeight {
  public:
      RunWeight(double alpha):sum(0.0),alpha(alpha){};
      inline void Update(double val) {
        sum=alpha*sum+val;
      }
    double sum;
  private:
    double alpha;
};

template <int tEMA=1,int tbias_corr=0>
class RunSum {
  public:
    RunSum(double alpha)
    :alpha(alpha),power_alpha(1.0),sum(0.0)
    {
    }
    inline void Update(double val)
    {
      if constexpr(tEMA)
        sum = alpha*sum + (1.0-alpha)*val;
      else
        sum = alpha*sum + val;

      if constexpr(tbias_corr)
        power_alpha*=alpha;
    }
    inline double Get() const
    {
      if constexpr(tbias_corr)
      {
        const double denom=1.0-power_alpha;
        assert(denom>0); // must call Update() before first Get()
        return sum/denom;
      } else
        return sum;
    }
  protected:
    double alpha,power_alpha,sum;
};

using RunSumEMA = RunSum<1,0>;
using RunSumEMA_BC = RunSum<1,1>;
using RunSumGEO = RunSum<0,0>;
using RunSumGEO_BC = RunSum<0,1>;

class RunMeanVar {
  public:
    RunMeanVar(double alpha)
    :alpha_(alpha),first_(true)
    {
      mean_ = var_ = 0.0;
    }
    void Update(double val)
    {
      /*if (first_) {
        mean_ = val;
        var_ = 0.0;
        first_ = false;
      } else*/
      {
        #if 0
          double delta = val - mean_;
          mean_ += (1 - alpha_) * delta;
          var_ = alpha_ * var_ + (1 - alpha_) * delta * delta;
        #else // slightly more stable
          // Welford
          double old_mean = mean_;
          mean_=alpha_*mean_+(1.0-alpha_)*val;
          var_=alpha_*var_+(1.0-alpha_)*((val-old_mean)*(val-mean_));
        #endif
      }
    }
    auto get()
    {
      return std::pair{mean_,std::max(0.0,var_)};
    }
  protected:
    double alpha_;
    bool first_;
    double mean_,var_;
};

namespace StrUtils {
      inline void StrUpper(std::string &str)
      {
        std::transform(str.begin(), str.end(),str.begin(), ::toupper);
      }
      inline std::string str_up(const std::string &str)
      {
        std::string ts=str;
        for (auto &c:ts) c=toupper(c);
        return ts;
      }
      inline void SplitToken(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters)
      {
        auto lastPos = str.find_first_not_of(delimiters, 0); // Skip delimiters at beginning.
        auto pos     = str.find_first_of(delimiters, lastPos); // Find first "non-delimiter".

        while (std::string::npos != pos || std::string::npos != lastPos)  {
          tokens.push_back(str.substr(lastPos, pos - lastPos)); // Found a token, add it to the vector.
          lastPos = str.find_first_not_of(delimiters, pos); // Skip delimiters.  Note the "not_of"
          pos = str.find_first_of(delimiters, lastPos); // Find next "non-delimiter"
        }
      }
      inline void RemoveWhite(std::string &str,const std::string &whites)
      {
        auto firstPos = str.find_first_not_of(whites);

        if (firstPos!=std::string::npos) {
          auto lastPos = str.find_last_not_of(whites);
          str=str.substr(firstPos,lastPos-firstPos+1);
        } else str="";
      }
      inline void SplitFloat(const std::string &str,std::vector<float>&x)
      {
        std::vector <std::string> tokens;
        SplitToken(str,tokens,",");
        for (auto &token:tokens) x.push_back(std::stof(token));
      }
};

namespace MathUtils {

  inline double calc_loglik_L1(double abs_e,double b)
  {
    return -std::log(2*b) - abs_e / b;
  }
  inline double calc_loglik_L2(double sq_e,double sigma2)
  {
    return -0.5*std::log(2*M_PI*sigma2) - 0.5 * sq_e / sigma2;
  }

  // inverse of pos. def. symmetric matrix
  class InverseSym
  {
    public:
      InverseSym(int n)
      :chol(n),n(n),b(n)
      {

      }
      void Solve(const vec2D &matrix,vec2D &sol,const double nu=0.0)
      {
        if (!chol.Factor(matrix,nu)) {
          for (int i=0;i<n;i++) {
             std::fill(std::begin(b),std::end(b),0.0);
             b[i]=1.0;
             chol.Solve(b,sol[i]);
          }
        };
      }
    protected:
      slmath::Cholesky chol;
      int n;
      vec1D b;
  };

  // estimate running covariance of vectors of len n
  class EstCov {
    public:
      EstCov(int n,double alpha=0.998,double init_val=1.0)
      :mcov(n,vec1D(n)),
      n(n),alpha(alpha)
      {
        for (int i=0;i<n;i++)
          mcov[i][i]=init_val;
      }
      void Update(const vec1D &x)
      {
        for (int j=0;j<n;j++) {
          for (int i=0;i<n;i++) {
            mcov[j][i] = alpha*mcov[j][i] + (1.0-alpha)*x[i]*x[j];
          }
        }
      }
      vec2D mcov;
      int n;
      double alpha;
  };


  template <typename T>
  T med3(T a,T b,T c)
  {
    if ((a<b && b<c) || (c<b && b<a)) {
      return b;
    } else if ((b < a && a < c) || (c < a && a < b)) {
      return a;
    } else
      return c;
  }

  inline int iLog2(int val) {
    int nbits=0;
    while (val>>=1) nbits++;
    return nbits;
  }
  inline double SumDiff(const std::vector<double> &v1,const std::vector<double> &v2)
  {
     if (v1.size()!=v2.size()) return -1;
     else {
       double sum=0.;
       for (size_t i=0;i<v1.size();i++) sum+=fabs(v1[i]-v2[i]);
       return sum;
     }
  }
  inline int32_t S2U(int32_t val)
  {
    if (val<0) val=2*(-val);
    else if (val>0) val=(2*val)-1;
    return val;
  }
  inline int32_t U2S(int32_t val)
  {
    if (val&1) val=((val+1)>>1);
    else val=-(val>>1);
    return val;
  }
  inline double norm2(const std::vector<double> &vec1,const std::vector<double> &vec2)
  {
     if (vec1.size()!=vec2.size()) return 0;
     else {
       double sum=0.;
       for (size_t i=0;i<vec1.size();i++) {double t=vec1[i]-vec2[i];sum+=t*t;};
       return sqrt(sum);
     }
  }
  inline double mean(const std::vector<double> &vec)
  {
    if (vec.size()) {
      double sum=std::accumulate(begin(vec),end(vec),0.0);
      return sum/static_cast<double>(vec.size());
    }
    return 0;
  }
  //contraharmonic mean
  inline double meanL(const std::vector<double> &vec)
  {
    if (vec.size()) {
      double sum0=0.0;
      double sum1=0.0;
      for (size_t i=0;i<vec.size();++i) {
        sum0+=(vec[i]*vec[i]);
        sum1+=vec[i];
      }
      if (sum1>0.0) return sum0 / sum1;
      else return 0.;
    }
    return 0.;
  }
  inline double linear_map_n(int n0,int n1,double y0,double y1,int idx)
  {
    double dx = static_cast<double>(n1-n0);
    double dy = y1-y0;
    return idx*(dy/dx)+y0;
  }
  template <typename T>
  T sgn(T x) {
    return (x > 0) - (x < 0);
    /*if (x>0) return 1;
    if (x<0) return -1;
    return 0;*/
  }
};

namespace miscUtils {

  enum class MapMode {rec,exp,tanh,power,sigmoid};

  template <MapMode mode>
  double decay_map(double gamma, double val)
  {
    if constexpr (mode == MapMode::rec) return 1.0 / (1.0 + gamma * val);
    else if constexpr (mode == MapMode::exp) return std::exp(-gamma * val);
    else if constexpr (mode == MapMode::tanh) return 1.0-std::tanh(gamma * val);
    else if constexpr (mode == MapMode::power) return std::pow(gamma, val);
    else if constexpr (mode == MapMode::sigmoid) return 1.0 / (1.0 + std::exp(gamma*(val-1.0)));
    return 0;
  }

/*static float rsqrt(float __x)
{
    float reciprocal;
    __asm__ __volatile__ (
        "movss %1, %%xmm0\n"
        "rsqrtss %%xmm0, %%xmm1\n"
        "movss %%xmm1, %0\n"
        :"=m"(reciprocal)
        :"m"(__x)
        :"xmm0", "xmm1"
    );
  return reciprocal;
}*/
  template <typename T>
  void swap_erase(std::vector<T>& e, std::size_t idx)
  {
    if (idx < e.size()) {
      std::swap(e[idx], e.back());
      e.pop_back();
    }
  }

  inline void RollBack(vec1D &data,double input)
  {
    if (data.size()) {
      std::memmove(&data[1],&data[0],(data.size()-1)*sizeof(double));
      data[0]=input;
    }
  }
  inline std::string getTimeStrFromSamples(int64_t numsamples,int64_t samplerate)
  {
   std::ostringstream ss;
   int h,m,s,ms;
   h=m=s=ms=0;
   if (numsamples>0 && samplerate>0) {
     while (numsamples >= 3600*samplerate) {++h;numsamples-=3600*samplerate;};
     while (numsamples >= 60*samplerate) {++m;numsamples-=60*samplerate;};
     while (numsamples >= samplerate) {++s;numsamples-=samplerate;};
     ms=round((numsamples*1000.)/samplerate);
   }
   ss << std::setfill('0') << std::setw(2) << h << ":" << std::setw(2) << m << ":" << std::setw(2) << s << "." << ms;
   return ss.str();
 }
 inline std::string getTimeStrFromSeconds(int seconds)
 {
   std::ostringstream ss;
   int h,m,s;
   h=m=s=0;
   if (seconds>0) {
      while (seconds >= 3600) {++h;seconds-=3600;};
      while (seconds >= 60) {++m;seconds-=60;};
      s=seconds;
   }
   ss << std::setfill('0') << std::setw(2) << h << ":" << std::setw(2) << m << ":" << std::setw(2) << s;
   return ss.str();
 }
 inline std::string ConvertFixed(double val,int digits)
 {
   std::ostringstream ss;
   ss << std::fixed << std::setprecision(digits) << val;
   return ss.str();
 }
};

namespace BitUtils {
  uint32_t get32HL(const uint8_t *buf);
  uint32_t get32LH(const uint8_t *buf);
  uint16_t get16LH(const uint8_t *buf);
  void put16LH(uint8_t *buf,uint16_t val);
  void put32LH(uint8_t *buf,uint32_t val);
  std::string U322Str(uint32_t val);

  inline int32_t count_bits32(uint32_t m)
  {
    #ifdef __GNUC__
      return m == 0 ? 0 : (32 - __builtin_clz(m));
    #else
      return std::bit_width(m);
    #endif
  }
}

#endif // UTILS_H
