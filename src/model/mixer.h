#ifndef MIXER_H
#define MIXER_H

#include "../global.h"
#include "model.h"
#include "domain.h"


class LogMixer
{
  public:
    enum class InitType {Zero,Uniform};
    LogMixer(int n,InitType type=InitType::Uniform);
    int Predict(const std::vector<int>&p);
    void Update(int bit,int rate);
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
