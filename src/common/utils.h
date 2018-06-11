#ifndef UTILS_H
#define UTILS_H

#include "../global.h"

// running exponential smoothing
// sum=alpha*sum+(1.0-alpha)*val, where 1/(1-alpha) is the mean number of samples considered
class RunExp {
  public:
      RunExp(double alpha):alpha(alpha){sum=0;};
      RunExp(double alpha,double sum):alpha(alpha),sum(sum){};
      inline double Get(){return sum;};
      inline void Update(double val) {
        sum=alpha*sum+(1.-alpha)*val;
      }
  private:
    double alpha,sum;
};

// running weighted sum: sum_{i=0}^n alpha^(n-i) val
class RunWeight {
  public:
      RunWeight(double alpha):alpha(alpha){sum=0;};
      inline double Get(){return sum;};
      inline void Update(double val) {
        sum=alpha*sum+val;
      }
  private:
    double alpha,sum;
};

class strUtils {
  public:
      static void strUpper(std::string &str)
      {
        std::transform(str.begin(), str.end(),str.begin(), ::toupper);
      }
};

class MathUtils {
  public:
      static int iLog2(int val) {
        int nbits=0;
        while (val>>=1) nbits++;
        return nbits;
      }
      static double SumDiff(const std::vector<double> &v1,const std::vector<double> &v2)
      {
         if (v1.size()!=v2.size()) return -1;
         else {
           double sum=0.;
           for (size_t i=0;i<v1.size();i++) sum+=fabs(v1[i]-v2[i]);
           return sum;
         }
      }
      static int32_t S2U(int32_t val)
      {
        if (val<0) val=2*(-val);
        else if (val>0) val=(2*val)-1;
        return val;
      }
      static int32_t U2S(int32_t val)
      {
        if (val&1) val=((val+1)>>1);
        else val=-(val>>1);
        return val;
      }
      static double L2Dist(const std::vector<double> &vec1,const std::vector<double> &vec2)
      {
         if (vec1.size()!=vec2.size()) return -1;
         else {
           double sum=0.;
           for (size_t i=0;i<vec1.size();i++) {double t=vec1[i]-vec2[i];sum+=t*t;};
           return sqrt(sum);
         }
      }
      static double L2Dist2(const std::vector<double> &vec1,const std::vector<double> &vec2)
      {
         if (vec1.size()!=vec2.size()) return -1;
         else {
           double sum=0.;
           for (size_t i=0;i<vec1.size();i++) {double t=vec1[i]-vec2[i];sum+=t*t;};
           return sum;
         }
      }
};

class miscUtils {
  public:
  // retrieve time string
  static std::string getTimeStrFromSamples(int numsamples,int samplerate)
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
 static std::string getTimeStrFromSeconds(int seconds)
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
 static std::string ConvertFixed(double val,int digits)
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
}

#endif // UTILS_H
