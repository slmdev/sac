#ifndef MIXER_H
#define MIXER_H

#include "../global.h"
#include "model.h"
#include "domain.h"

class LogMixer
{
  public:
    LogMixer(int n);
    int Predict(const std::vector <int>&p);
    void Update(int bit,int rate);
  private:
    inline int idiv_signed(int64_t val,int64_t s)
    {
      const int64_t half=int64_t{1} << (s-1);
      return val<0
        ? -(((-val)+ half)>>s)
        :   ((val  + half)>>s);
    };
    std::vector <int16_t>x;
    std::vector <int>w;
    int16_t pd;
    uint8_t n;
};

#endif
