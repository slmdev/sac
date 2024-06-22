#ifndef BIAS_H
#define BIAS_H

#include "../global.h"
#include "../common/utils.h"
#include "lms.h"

/*
gives a tiny gain
*/

#define h1y(v,k) (((v)>>k)^(v))

#define BIAS_ROUND_PRED 1

class BiasEstimator {
  class CntAvg {
    struct bias_cnt {
      int cnt;
      double val;
    };
    public:
      CntAvg(int isize,int scale=32,int freq0=4)
      :scale_(scale),freq0(freq0),bias(isize)
      {
        for (auto &x:bias) {x.cnt=0;x.val=0.0;};
      }
      double GetBias(int ctx)
      {
        double x=bias[ctx].val/double(bias[ctx].cnt+freq0);
        return x;
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
      int scale_,freq0;
      std::vector <bias_cnt>bias;
  };

  class CntExp {
    public:
      CntExp(int isize,double alpha)
      :bias(isize), alpha_(alpha)
      {
      }
      double GetBias(int ctx)
      {
        return bias[ctx];
      }
      void UpdateBias(int ctx, double delta) {
        bias[ctx] = alpha_ * bias[ctx] + (1.0-alpha_)*delta;
      }
    private:
      vec1D bias;
      double alpha_;
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
    BiasEstimator(double lms_mu=0.003,int cnt_scale=32,double nd_sigma=1.5,double nd_lambda=0.998)
    :
    mix_ada(32,SSLMS(4,lms_mu)),
    hist_input(8),hist_delta(8),
    bias(1<<20),
    cnt_freq(1<<20,cnt_scale), //,f0(1,cnt_scale),
    sigma(nd_sigma),lambda(nd_lambda)
    //sp(2,0.001)
    //nn(2,1)
    {
      ctx0=ctx1=ctx2=mix_ctx=0;
      p=0.0;
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

      //int c0=fabs(hist_delta[0])>t?1:0;
      //int c1=fabs(hist_delta[1])>t?1:0;
      //int c2=fabs(hist_delta[2])>t?1:0;

      ctx0=0;
      ctx0+=b0<<0;
      //ctx0+=b1<<1;
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

      const double bias0_a=cnt_freq.GetBias(ctx0);
      const double bias1_a=cnt_freq.GetBias(ctx1);
      const double bias2_a=cnt_freq.GetBias(ctx2);
      //double bias_mean = (bias0_a+bias1_a+bias2_a)/3.0;
      //double bias_med = med3(bias0_a,bias1_a,bias2_a);

      const double pbias=mix_ada[mix_ctx].Predict({bias0_a,bias1_a,bias2_a});
      //if (std::isnan(pbias)) std::cout << "nan";
      return pred+pbias;
    }
    void Update(double val) {
      #ifdef BIAS_ROUND_PRED
        const double delta=val-std::round(p);
      #else
        const double delta=val-p;
      #endif
      miscUtils::RollBack(hist_input,val);
      miscUtils::RollBack(hist_delta,delta);

      const double q=sigma*sqrt(var_est);
      const double lb=mean_est-q;
      const double ub=mean_est+q;

      if ( (delta>lb) && (delta<ub)) {
        cnt_freq.UpdateBias(ctx0,delta);
        cnt_freq.UpdateBias(ctx1,delta);
        cnt_freq.UpdateBias(ctx2,delta);
      }

      mix_ada[mix_ctx].Update(delta);

      mean_est=lambda*mean_est+(1.0-lambda)*delta;
      var_est=lambda*var_est+(1.0-lambda)*((delta-mean_est)*(delta-mean_est));

      //sp.Update(delta);
      //nn.Update(val);
    }
  private:
    std::vector<SSLMS> mix_ada;
    vec1D hist_input,hist_delta;
    int ctx0,ctx1,ctx2,mix_ctx;
    double p;
    //double alpha,p,bias0,bias1,bias2;
    vec1D bias;
    CntAvg cnt_freq;
    double mean_est,var_est,sigma,lambda;
};


#endif // BIAS_H
