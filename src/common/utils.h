#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
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

class StrUtils {
  public:
      static void StrUpper(std::string &str)
      {
        std::transform(str.begin(), str.end(),str.begin(), ::toupper);
      }
      static void SplitToken(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters)
      {
        auto lastPos = str.find_first_not_of(delimiters, 0); // Skip delimiters at beginning.
        auto pos     = str.find_first_of(delimiters, lastPos); // Find first "non-delimiter".

        while (std::string::npos != pos || std::string::npos != lastPos)  {
          tokens.push_back(str.substr(lastPos, pos - lastPos)); // Found a token, add it to the vector.
          lastPos = str.find_first_not_of(delimiters, pos); // Skip delimiters.  Note the "not_of"
          pos = str.find_first_of(delimiters, lastPos); // Find next "non-delimiter"
        }
      }
      static void RemoveWhite(std::string &str,const std::string &whites)
      {
        auto firstPos = str.find_first_not_of(whites);

        if (firstPos!=std::string::npos) {
          auto lastPos = str.find_last_not_of(whites);
          str=str.substr(firstPos,lastPos-firstPos+1);
        } else str="";
      }
      static void SplitFloat(const std::string &str,std::vector<float>&x)
      {
        std::vector <std::string> tokens;
        SplitToken(str,tokens,",");
        for (auto &token:tokens) x.push_back(std::stof(token));
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
      /*static double L2Dist2(const std::vector<double> &vec1,const std::vector<double> &vec2)
      {
         if (vec1.size()!=vec2.size()) return -1;
         else {
           double sum=0.;
           for (size_t i=0;i<vec1.size();i++) {double t=vec1[i]-vec2[i];sum+=t*t;};
           return sum;
         }
      }*/
};

class miscUtils {
  public:
  // retrieve time string
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
  static void RollBack(vec1D &data,double input)
  {
    if (data.size()) {
      for (int i=(int)(data.size()-1);i>0;i--)
        data[i]=data[i-1];
      data[0]=input;
    }
  }
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
