#ifndef MIXER_H
#define MIXER_H

#include <cmath>
#include "../global.h"
#include "model.h"
#include "domain.h"

template <int ExtraBits>
class LogMixerBase {
  public:
    LogMixerBase(int n,float minit=1.0)
    :n(n),x(n),w(n),pd(0)
    {
      if (n && minit>0.0)
      for (auto &weight:w)
        weight=static_cast<int>(std::lround(minit*wscale/n));
    }

    void Update(int bit,int rate)
    {
      int err=(bit<<PBITS)-pd;
      for (int i=0;i<n;i++)
      {
        int dw=idiv_signed(int64_t(x[i])*err*rate,ushift);
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
    static constexpr int wbits  = BM::WBITS + ExtraBits;
    static constexpr int wscale = 1 << wbits;
    static constexpr int wrange = 1 << (wbits + 3);
    static constexpr int ushift = BM::USHIFT - ExtraBits;

    static_assert(ExtraBits >= 0);
    static_assert(wbits < 28 && ushift > 0);

    inline int dot(std::span<int16_t>x,std::span<int>w)
    {
      int64_t sum=0;
      for (int i=0;i<n;i++)
        sum+=int64_t(w[i])*x[i];
      return idiv_signed(sum,wbits);
    }
    static inline int idiv_signed(int64_t val,int64_t s)
    {
      const int64_t half=int64_t{1} << (s-1);
      return val<0
        ? -(((-val)+ half)>>s)
        :   ((val  + half)>>s);
    };
    uint8_t n;
    std::vector <int16_t>x;
    std::vector <int>w;
    int16_t pd;
};

template <int ExtraBits=0>
class LogMixerExtra : public LogMixerBase<ExtraBits>
{
  using Base=LogMixerBase<ExtraBits>;
  using Base::n,Base::dot,Base::x,Base::w,Base::pd;
  public:
    LogMixerExtra(int n,float minit=1.0)
    :Base(n,minit){};
    int Predict(std::span<const int>p)
    {
      for (int i=0;i<n;i++)
        x[i]=myDomain.Fwd(p[i]);
      const int logit=dot(x,w);
      return (pd=myDomain.Inv(logit));
    }
    int Predict(std::initializer_list<int>p)
    {
      return Predict(std::span<const int>(p.begin(), p.size()));
    }
};

using LogMixer=LogMixerExtra<0>;

class LogMixerBias : public LogMixerBase<5>
{
  public:
    LogMixerBias(int n,float minit=0.0)
    :LogMixerBase(n,minit){};
    int Predict(int p1,std::span<const int>p);
    int Predict(int p1,std::initializer_list<int>p)
    {
      return Predict(p1,std::span<const int>(p.begin(), p.size()));
    }
};

#endif
