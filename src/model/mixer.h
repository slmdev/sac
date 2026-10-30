#ifndef MIXER_H
#define MIXER_H

#include <cmath>
#include "../global.h"
#include "model.h"
#include "domain.h"


class LogMixer
{
  public:
    enum class InitType {Zero,Uniform};
    LogMixer(int n,InitType type=InitType::Uniform);
    int Predict(std::span<const int>p);
    int Predict(std::initializer_list<int>p)
    {
      return Predict(std::span<const int>(p.begin(), p.size()));
    }
    void Update(int bit,int rate);
    inline void SetConfidence(float s)
    {
      for (auto &weight:w)
        weight = static_cast<int>(std::lround(s*WSCALE));
    }
    inline int idiv_signed(int64_t val,int64_t s)
    {
      const int64_t half=int64_t{1} << (s-1);
      return val<0
        ? -(((-val)+ half)>>s)
        :   ((val  + half)>>s);
    };
    inline int dot(std::span<int16_t>x,std::span<int>w)
    {
      int64_t sum=0;
      for (int i=0;i<n;i++)
        sum+=int64_t(w[i])*x[i];
      return idiv_signed(sum,WBITS);
    }
    std::vector <int16_t>x;
    std::vector <int>w;
    int16_t pd;
    uint8_t n;
};

#endif
