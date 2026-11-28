#ifndef MIXER_H
#define MIXER_H

#include <cmath>
#include "../global.h"
#include "model.h"
#include "domain.h"

//x = input in log-domain
//w = weights (precision=BM::WBITS+ExtraBits
template <int ExtraBits>
class LogMixerBase {
  public:
    //minit = initial mixer confidence 0.0=none,1.0=uniform
    LogMixerBase(int n,float minit=1.0)
    :n(n),w(n)
    {
      if (n<=0)
        throw std::invalid_argument("LogMixerBase: invalid inputsize");
      if (minit>0.0)
        for (auto &weight:w)
          weight=static_cast<int>(std::lround(minit*wscale/n));
    }
    int Predict(std::span<const int16_t>x)
    {
      const int logit=dot(x,w);
      return (pd=myDomain.Inv(logit));
    }
    void Update(std::span<const int16_t> x,int bit,int rate)
    {
      const int err=(bit<<PBITS)-pd;
      const int gain=int64_t(err)*rate;
      for (int i=0;i<n;i++)
      {
        int dw=idiv64(int64_t(x[i])*gain,ushift);
        w[i]=std::clamp(w[i]+dw,-wrange,wrange-1);
      }
    }
    void SetConfidence(float s)
    {
      for (int i=0;i<n;i++)
        SetConfidence(i,s/n);
    }
    void SetConfidence(int idx,float s)
    {
      w[idx]=static_cast<int>(std::lround(s*wscale));
    }

  protected:
    static int dot(std::span<const int16_t>x,std::span<const int>w)
    {
      int64_t sum=0;
      for (std::size_t i=0;i<x.size();++i)
        sum+=int64_t(w[i])*x[i];
      return idiv64(sum,wbits);
    }
    int16_t n;
    std::vector <int>w;
    int16_t pd=0;

  private:
    static constexpr bool branchless_div=true;
    static constexpr int wbits  = BM::WBITS + ExtraBits;
    static constexpr int wscale = 1 << wbits;
    static constexpr int wrange = 1 << (wbits + 3);
    static constexpr int ushift = BM::USHIFT - ExtraBits;

    static_assert(ExtraBits >= 0);
    static_assert(wbits < 28 && ushift > 0);

    static int idiv64(int64_t val,int s)
    {
      const int64_t half=int64_t{1} << (s-1);
      if constexpr(branchless_div) {
        return static_cast<int>((val+half+(val>>63))>>s);
      } else {
        return static_cast<int>(val<0
          ? -(((-val)+ half)>>s)
          :   ((val  + half)>>s));
      }
    }
};

template <int ExtraBits=0>
class LogMixerExtra : public LogMixerBase<ExtraBits>
{
  using Base=LogMixerBase<ExtraBits>;
  using Base::n,Base::dot,Base::w,Base::pd;
  public:
    LogMixerExtra(int n,float minit=1.0)
    :Base(n,minit),x(n){};
    int Predict(std::span<const int>p) {
      for (int i=0;i<n;i++)
        x[i]=myDomain.Fwd(p[i]);
      return Base::Predict(x);
    }
    int Predict(std::initializer_list<int>p) {
      return Predict(std::span<const int>(p.begin(), p.size()));
    }
    void Update(int bit,int rate) {
      Base::Update(x,bit,rate);
    }
  private:
    std::vector<int16_t> x;
};

class LogMixerBias : public LogMixerBase<5>
{
  using Base=LogMixerBase<5>;
  public:
    LogMixerBias(int n,float minit=0.0)
    :LogMixerBase(n,minit),x(n){};
    int Predict(int p1,std::span<const int>p);
    int Predict(int p1,std::initializer_list<int>p) {
      return Predict(p1,std::span<const int>(p.begin(), p.size()));
    }
    void Update(int bit,int rate) {
      Base::Update(x,bit,rate);
    }
  private:
    std::vector<int16_t> x;
};

using LogitMixer=LogMixerBase<0>;
using LogMixer=LogMixerExtra<0>;

#endif
