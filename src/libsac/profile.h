#ifndef PROFILE_H
#define PROFILE_H

#include <vector>
#include "map.h"

class SacProfile {
  public:
    struct FrameStats {
      int maxbpn,maxbpn_map;
      bool enc_mapped;
      int32_t blocksize,minval,maxval,mean;
      Remap mymap;
    };
    struct coef {
      float vmin,vmax,vdef;
    };
      SacProfile():type(0){};
      void Init(int numcoefs,int ptype)
      {
         coefs.resize(numcoefs);
         type=ptype;
      }
      SacProfile(int numcoefs,int type)
      :type(type),coefs(numcoefs)
      {

      }
      void Set(int num,double vmin,double vmax,double vdef)
      {
        if (num>=0 && num< static_cast<int>(coefs.size())) {
          coefs[num].vmin=vmin;
          coefs[num].vmax=vmax;
          coefs[num].vdef=vdef;
        }
      }
      void Set(int num,const std::vector<float>&x)
      {
        if (num>=0 && num<static_cast<int>(coefs.size()) && (x.size()>=3)) {
          coefs[num].vmin=x[0];
          coefs[num].vmax=x[1];
          coefs[num].vdef=x[2];
        }
      }
      double Get(int num) const {
        if (num>=0 && num< static_cast<int>(coefs.size())) {
            return coefs[num].vdef;
        } else return 0.;
      }
      int type;
      std::vector <coef> coefs;
};

#endif // PROFILE_H
