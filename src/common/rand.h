#ifndef RAND_H
#define RAND_H

#include "../global.h"

class Random {
  public:
    Random():engine(time(0)){};
    Random(uint32_t seed):engine(seed){};
    double r_01() { // [0,1)
      return std::uniform_real_distribution<double>{0,1}(engine);
    };
    double r_01open() { // [0,1)
      return std::uniform_real_distribution<double>{std::nextafter(0.0, std::numeric_limits<double>::max()),1.0}(engine);
    };
    double r_01closed() { // [0,1]
      return std::uniform_real_distribution<double>{0,std::nextafter(1.0, std::numeric_limits<double>::max())}(engine);
    };
    double r_int(double imin,double imax) { //double in [imin,imax]
      return std::uniform_real_distribution<double>{imin,std::nextafter(imax, std::numeric_limits<double>::max())}(engine);
    };
    uint32_t ru_int(uint32_t imin,uint32_t imax) { //int in [imin,imax]
      return std::uniform_int_distribution<uint32_t>{imin, imax}(engine);
    };
    double r_norm(double mu,double sigma) { // normal
      return std::normal_distribution<double>{mu,sigma}(engine);
    }
    double r_lognorm(double mu,double sigma) { // log-normal
      return exp(std::normal_distribution<double>{mu,sigma}(engine));
    }
    uint32_t ru_geo(double p) { // geometric
      return std::geometric_distribution<uint32_t>{p}(engine);
    }
    uint32_t ru_poi(double lambda) { // poisson
      return std::poisson_distribution<uint32_t>{lambda}(engine);
    }
    bool event(double p) {
      if (r_01()<p) return true;
      else return false;
    };
  private:
    std::mt19937 engine;
};

#endif
