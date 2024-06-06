#ifndef MIXER_H
#define MIXER_H

#include "model.h"
#include "domain.h"

// adaptive linear 2-input mix
// maximum weight precision 16-Bit
class Mix2
{
  public:
    virtual int Predict(int _p1,int _p2)=0;
    virtual void Update(int bit,int rate)=0;
};

class Mix2Linear : public Mix2
{
  public:
    Mix2Linear(){Init(WSCALEh);};
    void Init(int iw){w=iw;};
    //pm=(1-w)*p1+w*p2
    int Predict(int _p1,int _p2)
    {
      p1=_p1;p2=_p2;
      pm = p1+idiv_signed32((p2-p1)*w,WBITS);
      pm = clamp(pm,1,PSCALEm);
      return pm;
    }
    int w,p1,p2,pm;
  protected:
    inline int idiv_signed32(int val,int s){return val<0?-(((-val)+(1<<(s-1)))>>s):(val+(1<<(s-1)))>>s;};
    inline void upd_w(int d,int rate) {int wd=idiv_signed32(rate*d,PBITS);w=clamp(w+wd,0,int(WSCALE));};
};

class Mix2LeastSquares : public Mix2Linear  {
  public:
    // w_(i+1)=w_i + rate*(p2-p1)*e
    void Update(int bit,int rate)
    {
      int e=(bit<<PBITS)-pm;
      int d=idiv_signed32((p2-p1)*e,PBITS);
      upd_w(d,rate);
    }
};

class Mix2LeastCost : public Mix2Linear {
  public:
    void Update(int bit,int rate)
    {
      int d;
      //if (bit) d=(((p2-p1)<<PBITS)*uint64_t(div32tbl[pm]))>>32;
      //else d=(((p1-p2)<<PBITS)*uint64_t(div32tbl[PSCALE-pm]))>>32;
      if (bit) d=((p2-p1)<<PBITS)/pm;
      else d=((p1-p2)<<PBITS)/(PSCALE-pm);
      upd_w(d,rate);
    }
};

class NMixLogistic
{
  enum {WRANGE=1<<19};
  std::vector <int16_t> x;
  std::vector <int>w;

  int16_t pd;
  uint8_t n;
  public:
    NMixLogistic(int n)
    :x(n),w(n),pd(0),n(n)
     {
       Init(0);
     };
    void Init(int iw){
      for (int i=0;i<n;i++) w[i]=iw;
    };
    int Predict(const std::vector <int>&p)
    {
      int64_t sum=0;
      for (int i=0;i<n;i++)
      {
         x[i]=myDomain.Fwd(p[i]);
         sum+=int64_t(w[i]*x[i]);
      }
      sum=idiv_signed64(sum,WBITS);
      pd=clamp(myDomain.Inv(sum),1,PSCALEm);
      return pd;
    }
    void Update(int bit,int rate)
    {
      int err=(bit<<PBITS)-pd;
      for (int i=0;i<n;i++)
      {
         int de=idiv_signed32(x[i]*err,myDomain.dbits);
         upd_w(i,idiv_signed32(de*rate,myDomain.dbits));
      }
    };
  protected:
    inline int idiv_signed32(int val,int s){return val<0?-(((-val)+(1<<(s-1)))>>s):(val+(1<<(s-1)))>>s;};
    inline int idiv_signed64(int64_t val,int64_t s){return val<0?-(((-val)+(1<<(s-1)))>>s):(val+(1<<(s-1)))>>s;};
    inline void upd_w(int i,int wd){w[i]=clamp(w[i]+wd,-WRANGE,WRANGE-1);}
};

#endif // MIXER_H
