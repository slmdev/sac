#ifndef PMODEL_H
#define PMODEL_H

#include <cmath>
#include "../global.h"
#include "../model/model.h"

using namespace BM;

//multiplicative correction model for standard laplace
class ScaleModel {
  static constexpr double MAX_LOG_STEP=std::log(2.0);
  static constexpr double MAX_BIAS=std::log(2.0);
  static constexpr double MAX_CORR=20.0;
  static constexpr double MEAN_EPS=1E-3;
public:
  ScaleModel(double scale_rate=0.003,double bias_rate=1E-4)
  :scale_rate(scale_rate),bias_rate(bias_rate)
  {}
  void Begin(double corr_retention=1.0,double bias_retention=0.5)
  {
    if (corr_retention<1.0)
      corr *= corr_retention;
    if (bias_retention<1.0)
      bias *= bias_retention;
  }
  int Predict(double mean, uint32_t plane)
  {
    const double c=std::exp(corr);
    //r = std::pow(plane / (mean*c+MEAN_EPS),0.75);
    r = static_cast<double>(plane)/(mean*c+MEAN_EPS);
    const double t = std::exp(-r+bias);
    ps = t / (1.0 + t);
    return std::clamp(
      static_cast<int>(std::lround(ps*PSCALE)),1,PSCALEm);
  }
  void Update(int bit)
  {
    const double error = (bit-ps);
    if (scale_rate>0.0) corr = std::clamp(corr+std::clamp(scale_rate*error*r,-MAX_LOG_STEP,MAX_LOG_STEP),-MAX_CORR,MAX_CORR);
    if (bias_rate>0.0) bias = std::clamp(bias+bias_rate*error,-MAX_BIAS,MAX_BIAS);
  }
private:
  double r=0.0,ps=0.0;
  double corr=0.0,bias=0.0;
  double scale_rate,bias_rate;
};

class LaplaceModel {
  public:
    static int PredictLogit(double mean,uint32_t plane)
    {
      const double z=-static_cast<double>(plane)/(mean+1E-5);
      const int zlog=std::lround(std::clamp(z*LSCALE,
      static_cast<double>(BM::DMIN),static_cast<double>(BM::DMAX)));
      return zlog;
    }
};

#endif

