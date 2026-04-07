#ifndef HM_H
#define HM_H

#include <cmath>
#include <vector>
#include <type_traits>
#include <algorithm>
#include "ls.h"

/*
  Hierarchical mixer using different gatings f
  standard LMS = HMix<L2,Gate0>

  Gate0: y := y + b*x
  Gate1: y := y + f(a)*b*x //tanh
  Gate2: y := y + f(a*ctx)*b*x
  Gate3: y := (1.0+scale*f(a))*y + b*x //tanh with small scale
*/

namespace HM {

//Activation functions
template<typename T>
concept Activation = requires(double x) {
  { T::f(x) } -> std::same_as<double>;
  { T::df(x) } -> std::same_as<double>;
};
struct Sigmoid { //[0,1]
  static double f(double x){return 1.0/(1.0+std::exp(-x));}
  static double df(double x){double s=f(x);return 1.0*s*(1.0-s);}
};
struct Tanh { //[-1,1]
  static double f(double x){return std::tanh(x);}
  static double df(double x){double t = std::tanh(x);return 1.0 - t * t;}
};
struct TanhG {//scaled and shifted tanh [0,1]
  static double f(double x) {return 0.5*(std::tanh(x)+1.0);}
  static double df(double x){return 0.5*(1.0 - std::tanh(x)*std::tanh(x));}
};
struct ReLU {
  static double f(double x){return x>0.0 ? x : 0.0;}
  static double df(double x){return x>0.0 ? 1.0 : 0.0;}
};
struct LeakyReLU {
  static constexpr double alpha = 0.001;
  static double f(double x){return x > 0.0 ? x : alpha * x;}
  static double df(double x){return x > 0.0 ? 1.0 : alpha;}
};


// Gate definitions
struct GateParam{double a,b;};
struct GateInput{double prev,x,ctx;};

struct Gate0 {
  static double fw(const GateInput &gi,const GateParam &gp) {
    return gi.prev + inject(gi,gp)*gi.x;
  }
  static GateParam grad(double delta,const GateInput &gi,const GateParam&)
  {
    return {0.0,delta*gi.x};
  }
  static constexpr double carry(double){return 1.0;}
  static double inject(const GateInput&,const GateParam &gp){return gp.b;}
};

template<Activation A>
struct Gate1 {
  static double fw(const GateInput &gi,const GateParam &gp) {
    return gi.prev + inject(gi,gp)*gi.x;
  }
  static GateParam grad(double delta,const GateInput &gi,const GateParam &gp)
  {
    double ga = delta * A::df(gp.a) * gp.b * gi.x;  // dL/da[i]
    double gb = delta * A::f(gp.a) * gi.x;     // dL/db[i]
    return {ga,gb};
  }
  static constexpr double carry(double){return 1.0;}
  static double inject(const GateInput&,const GateParam &gp){

    return A::f(gp.a) * gp.b;
  }
};

template<Activation A>
struct Gate2 {
  static double fw(const GateInput &gi,const GateParam &gp) {
    return gi.prev + inject(gi,gp)*gi.x;
  }
  static GateParam grad(double delta,const GateInput &gi,const GateParam &gp)
  {
    double z = std::clamp(gi.ctx*gp.a,-5.0,5.0);
    //double z = gi.ctx*gp.a;
    double ga = delta * gp.b * gi.x * gi.ctx * A::df(z);  // dL/da[i]
    double gb = delta * A::f(z) * gi.x;     // dL/db[i]
    return {ga,gb};
  }
  static constexpr double carry(double)  {return 1.0;}
  static double inject(const GateInput &gi,const GateParam &gp){
    //double z = gi.ctx*gp.a;
    double z = std::clamp(gp.a*gi.ctx,-5.0,5.0);
    return A::f(z)*gp.b;
  }
};

template<Activation A>
struct Gate3 {
  static constexpr double carry_scale=0.05; //5% around carry=1.0
  static double fw(const GateInput &gi,const GateParam &gp) {
    return carry(gp.a)*gi.prev + inject(gi,gp)*gi.x;
  }
  static GateParam grad(double delta,const GateInput &gi,const GateParam &gp)
  {
    double ga = delta * gi.prev * A::df(gp.a) * carry_scale; // dL/da[i]
    double gb = delta * gi.x;     // dL/db[i]
    return {ga,gb};
  }
  static double carry(double a) {
    return 1.0 + carry_scale * A::f(a);
  };
  static double inject(const GateInput&,const GateParam &gp) {return gp.b;}
};

//Regularization
template<typename T>
concept WeightDecay = requires(double w) {
  { T::apply(w) } -> std::same_as<double>;
};
struct NoReg {static double apply(double) {return 0.0;}};
template <double nu=0.001>
struct L1Reg {static double apply(double w) {
  return nu*MathUtils::sgn(w);}
};
template <double nu=0.001>
struct L2Reg {static double apply(double w) {
  return nu*(w);}
};

//detect gate2
template<typename T>
struct is_gate2 : std::false_type {};

template<Activation A>
struct is_gate2<Gate2<A>> : std::true_type {};

template<Loss::LossFunction LF,class Gate,WeightDecay Reg=NoReg,int init_type=0>
class HMix : public LS {
public:
  HMix(int n,double mu,double beta=0.95,double gate_init=0.0)
  :LS(n,mu),
  input(n, {.prev = 0.0, .x = 0.0, .ctx = 0.0}),
  param(n, {.a = gate_init, .b = 0.0}),
  stats(n, {.a = 0, .b = 0}),
  rctx(n,0.95),beta(beta)
  {
    double init_b=0.0;
    if constexpr(init_type==1)
      init_b = 1.0;
    else if constexpr(init_type==2)
      init_b = 1.0/n;

    for (int i=0;i<n;++i) {
      if constexpr(init_type==3)
        init_b = 1.0/(1.0+i);
      param[i].b = init_b;
    }
  }

  double Predict(span_cf64 x) override
  {
    double pred=0.0;
    for (std::size_t i = 0;i<n;++i) {
      double mctx=0.0;

      if constexpr (is_gate2<Gate>::value)
      {
        double tctx=std::log1p(std::abs(x[i]));
        rctx[i].Update(tctx);
        auto [m,v] = rctx[i].Get();
        //mctx = m/(std::sqrt(v)+1E-5);
        mctx = m;
      }

      input[i] = {.prev=pred,.x=x[i],.ctx=mctx}; //set gate input
      pred = Gate::fw(input[i],param[i]);
    }
    return pred;
  }

  void Update(span_cf64,double error) override
  {
    // dL/d(out[n-1])
    double delta=LF::grad(error);
    //backpropagation
    for (int i=n - 1;i >= 0;--i) {
      // use pre-update value for backprop
      const double a_old = param[i].a;

      auto g = Gate::grad(delta,input[i],param[i]);

      stats[i].a=beta*stats[i].a+(1.0-beta)*(g.a*g.a);
      stats[i].b=beta*stats[i].b+(1.0-beta)*(g.b*g.b);
      double sga = g.a/(std::sqrt(stats[i].a)+SACCfg::LMS_ADA_EPS);
      double sgb = g.b/(std::sqrt(stats[i].b)+SACCfg::LMS_ADA_EPS);

      // update weights including regularization
      param[i].a += mu*(sga);;
      param[i].b += mu*(sgb - Reg::apply(param[i].b));

      // dL/d(out[i-1]) = dL/d(out[i]) * carry
      delta *= Gate::carry(a_old);
    }
  }

  // get the weighted contribution of stage i on the final prediction
  double GetWeight(int i) const override
  {
    double w = Gate::inject(input[i],param[i]);
    for (std::size_t j=i+1;j<n;++j)
      w *= Gate::carry(param[j].a);
    return w;//std::max(w,0.0);
  }

private:
  std::vector<GateInput> input;
  std::vector<GateParam> param;
  std::vector<GateParam> stats;
  std::vector<RunMeanVar> rctx;

  double beta;
};

};

#endif

