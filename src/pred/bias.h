#ifndef BIAS_H
#define BIAS_H

#include "../global.h"
#include "../common/utils.h"
#include "lms.h"

/*
gives a tiny gain
*/

class BiasEstimator {
  class BiasCnt {
    struct bias_cnt {
      int cnt;
      double val;
    };
    public:
      BiasCnt(int scale,int isize)
      :scale_(scale),bias(isize)
      {
        for (auto &x:bias) {x.cnt=0;x.val=0.0;};
      }
      double GetBias(int ctx)
      {
        return bias[ctx].val/(bias[ctx].cnt+1);
      }
      void UpdateBias(int ctx,double delta) {
        bias[ctx].val+=delta;
        bias[ctx].cnt++;
        if (bias[ctx].cnt>scale_) {
          bias[ctx].val/=2.0;
          bias[ctx].cnt/=2;
        }
      }
    private:
      int scale_;
      std::vector <bias_cnt>bias;
  };

  double med3(double a,double b, double c) {
    if ((a<b && b<c) || (c<b && b<a)) {
      return b;
    } else if ((b < a && a < c) || (c < a && a < b)) {
      return a;
    } else
      return c;
  }

  public:
    BiasEstimator(double mu=0.002,double lambda=0.998,int re_scale=32)
    :mix_ada(32,SSLMS(3,mu)),
    hist_input(8),hist_delta(8),
    bias(1<<20),
    Bias0(re_scale,1<<20),lambda(lambda)
    {
      ctx0=ctx1=ctx2=mix_ctx=0;
      p=0.0;
      //lambda=0.998;
      mean_est=var_est=0.;
    }
    void CalcContext()
    {
      int b0=hist_input[0]>p?0:1;
      //int b1=hist_input[1]>p?0:1;
      int b2=hist_delta[0]<0?0:1;
      int b3=hist_delta[1]<0?0:1;
      int b4=hist_delta[2]<0?0:1;
      //int b42=hist_delta[3]<0?0:1;
      int b5=hist_delta[1]<hist_delta[0]?0:1;
      int b6=hist_delta[2]<hist_delta[1]?0:1;
      int b7=hist_delta[3]<hist_delta[2]?0:1;
      int b8=hist_delta[4]<hist_delta[3]?0:1;
      int b9=(fabs(hist_delta[0]))>32?0:1;
      int b10=2*hist_input[0]-hist_input[1]>p?0:1;
      int b11=3*hist_input[0]-3*hist_input[1]+hist_input[2]>p?0:1;

      double t=(fabs(hist_delta[0])+fabs(hist_delta[1])+fabs(hist_delta[2])+fabs(hist_delta[3])+fabs(hist_delta[4]))/5.;
      mix_ctx=0;
      if (t>512) mix_ctx=2;
      else if (t>32) mix_ctx=1;
      else mix_ctx=0;

      //int c0=fabs(hist_delta[0])>t?1:0;
      //int c1=fabs(hist_delta[1])>t?1:0;
      //int c2=fabs(hist_delta[2])>t?1:0;

      ctx0=0;
      /*ctx0+=c0<<0;
      ctx0+=c1<<1;
      ctx0+=c2<<2;*/
      ctx0+=b0<<0;
      ctx0+=b2<<1;
      ctx0+=b9<<2;
      ctx0+=b10<<3;
      ctx0+=b11<<4;

      ctx1=64;
      ctx1+=b2<<0;
      ctx1+=b3<<1;
      ctx1+=b4<<2;

      ctx2=1024;
      ctx2+=b5<<0;
      ctx2+=b6<<1;
      ctx2+=b7<<2;
      ctx2+=b8<<3;
      //ctx2+=mix_ctx<<4;
    }
    double Predict(double pred)
    {
      p=pred;
      CalcContext();

      //bias0=bias[ctx0];
      //bias1=bias[ctx1];
      //bias2=bias[ctx2];

      double bias0_a=Bias0.GetBias(ctx0);
      double bias1_a=Bias0.GetBias(ctx1);
      double bias2_a=Bias0.GetBias(ctx2);
      //double bias_mean = (bias0_a+bias1_a+bias2_a)/3.0;
      //double bias_med = med3(bias0_a,bias1_a,bias2_a);

      double pbias=mix_ada[mix_ctx].Predict({bias0_a,bias1_a,bias2_a});
      //if (std::isnan(pbias)) std::cout << "nan";

      return pred+pbias;
    }
    void Update(double val) {
      const double delta=val-p;
      miscUtils::RollBack(hist_input,val);
      miscUtils::RollBack(hist_delta,val-p);

      double sigma=1.5;
      double lb=mean_est-sigma*sqrt(var_est);
      double ub=mean_est+sigma*sqrt(var_est);
      bool is_in=delta>lb && delta<ub;

      if (is_in) {
        /*bias[ctx0]=alpha*bias[ctx0]+(1.0-alpha)*delta;
        bias[ctx1]=alpha*bias[ctx1]+(1.0-alpha)*delta;
        bias[ctx2]=alpha*bias[ctx2]+(1.0-alpha)*delta;*/

        Bias0.UpdateBias(ctx0,delta);
        Bias0.UpdateBias(ctx1,delta);
        Bias0.UpdateBias(ctx2,delta);
      }


      mix_ada[mix_ctx].Update(delta);

      mean_est=lambda*mean_est+(1.0-lambda)*delta;
      var_est=lambda*var_est+(1.0-lambda)*((delta-mean_est)*(delta-mean_est));
    }
  private:
    std::vector<SSLMS> mix_ada;
    vec1D hist_input,hist_delta;
    int ctx0,ctx1,ctx2,mix_ctx;
    double p;
    //double alpha,p,bias0,bias1,bias2;
    vec1D bias;
    BiasCnt Bias0;
    double mean_est,var_est,lambda;
};


#endif // BIAS_H
