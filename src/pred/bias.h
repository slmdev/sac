#ifndef BIAS_H
#define BIAS_H

#include "../global.h"

class BC {
  public:
      BC():hist(10),cmix(2,0.001)
      {
        for (int i=0;i<4;i++) {esum[i]=0;en[i]=0;};
      }
    double Predict(double p) {
      ctx=0;
      if (hist[0]>p) ctx+=1;
      if (hist[1]>p) ctx+=2;
      //ctx=0;

      if (en[ctx]) {
        bias=esum[ctx]/double(en[ctx]);
        //cout << bias << " ";
      } else bias=0;
      pm=p+bias;

      std::vector<double> pred={p,pm};
      double px=cmix.Predict(pred);
      //cout << pm << " " << px << endl;
      return px;
    }
    void Update(double val)
    {
       cmix.Update(val);
       double e=val-pm;
       double ae=fabs(e);
       if (ae<256) {
         esum[ctx]+=e;
         en[ctx]++;
         if (en[ctx]>32) {
            en[ctx]=en[ctx]/2;
            esum[ctx]=esum[ctx]/2.;
         }
       }
       for (int i=9;i>0;i--) hist[i]=hist[i-1];hist[0]=val;
    }
  private:
    double pm,bias;
    double esum[4];
    int en[4];
    std::vector<double> hist;
    int ctx;
    SSLMS cmix;
};

class BiasCorrection {
  public:
    struct ectx {
      int sum;
      int num;
      int succ;
    };
    BiasCorrection():hist(10),err(512) {
      ctx=0;lerr=0;
    };
    int32_t Predict(int pred)
    {
      ctx=0;
      if (hist[0]>pred) ctx+=1;
      if (hist[1]>pred) ctx+=2;
      //if ((2*hist[0]-hist[1])>pred) ctx+=4;
      //if ((2*hist[1]-hist[2])>pred) ctx+=8;
      //if (lerr>0) ctx+=16;

      p=pred;
      bias=0;

      if (err[ctx].num) {
         bias=(err[ctx].sum)/err[ctx].num;
         if (abs(bias)) {
           double r=err[ctx].succ/(double)err[ctx].num;
           //cout << "(" << r << "," << bias << ")";
           if (r>0.9) return (p+bias);
           else return p;
         } else return p;
      } else return p;
    }
    void Update(int32_t val) {
      int error=val-(p+bias);
      err[ctx].sum+=error;
      err[ctx].num++;
      if (abs(val-(p+bias))<=abs(val-p)) err[ctx].succ++;
      if (err[ctx].num==2048) {
        err[ctx].sum/=2;
        err[ctx].num/=2;
        err[ctx].succ/=2;
      }
      for (int i=9;i>0;i--) hist[i]=hist[i-1];hist[0]=val;
      lerr=error;
    }
  protected:
    std::vector<double> hist;
    std::vector <ectx>err;
    int ctx,p,bias,lerr;

};


#endif // BIAS_H
