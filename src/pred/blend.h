#ifndef BLEND_H
#define BLEND_H

//#include <memory>

// blend two expert outputs via sigmoid
class Blend2 {
  public:
    Blend2(double beta=0.95,double theta0=1.0,double theta1=0.0,double scale=5.0)
    :w(0.5),
     th0(theta0),th1(theta1),scale(scale),
     rsum(beta)
    {
    }
    double Predict(double p0,double p1) const
    {
      return w*p0 + (1.0-w)*p1;
    }

    void Update(double score0,double score1)
    {
      rsum.Update(score1-score0);

      double delta=rsum.Get();
      double z = th0*delta + th1;
      z = std::clamp(z,-scale,scale);
      w = 1.0 / (1.0 + std::exp(-z));
   }
  private:
    double w,th0,th1,scale;
    RunSumEMA rsum;
};


class BlendRegret
{
  static constexpr double EPS=1E-8;
  public:
    BlendRegret(std::size_t n,double alpha,double beta)
    :n_(n),beta(beta),px(0.0),
    x(n),w(n),zm(n),
    rsum(n,RunMeanVar(alpha))
    {
      if (n)
        std::fill(begin(w),end(w),1.0/n); //init equal weight

    };
    double Predict(const vec1D &input)
    {
      x=input;
      px=slmath::dot(x,w);
      return px;
    }
    void UpdateRegret(double target)
    {
      double loss_px = std::abs(target-px);
      for (std::size_t i=0;i<n_;i++) {
        double loss_pi=std::abs(target-x[i]);
        // if regret < 0 -> expert better then blend
        double regret=(loss_pi - loss_px);
        rsum[i].Update(regret);
      }
    }
    // softmax w_i = exp(-beta * normalized_regret)
    void UpdateWeights()
    {
      double max_z = -std::numeric_limits<double>::infinity();
      for (std::size_t i=0;i<n_;i++) {
        auto [mean,var] = rsum[i].Get(); //regret

        // scaled signal-to-noise
        zm[i]= -beta*mean/(std::sqrt(var)+EPS);
        max_z = std::max(max_z,zm[i]);
      }

      //best expert has highest z-score -> weight=exp(0)=1
      double total=0.0;
      for (std::size_t i=0;i<n_;i++) {
        w[i] = std::exp(zm[i]-max_z);
        total += w[i];
      }
      //normalize weights, total >= 1 from max-trick
      const double inv_total=1.0/total;
      for (double &val : w) val *= inv_total;
    }
    const vec1D &Weights()const {return w;}
  private:
    std::size_t n_;
    double beta,px;
    vec1D x,w,zm;
    std::vector <RunMeanVar> rsum;
};

      #if 0
        double hr=0.5,beta_nu=0.05;
        const double Hmax = std::log(np);
        const double Htarget = Hmax*hr;
        double H = 0.0;
        for (std::size_t i=0;i<np;i++) {
          if (w[i]>1E-8) H-=w[i]*log(w[i]);
        }
        beta += beta_nu*(H-Htarget);
        beta = std::clamp(beta,1.0,15.0);
        //std::cout << H << ' ' << Htarget << ' ' << beta << '\n';
      #endif


#endif

