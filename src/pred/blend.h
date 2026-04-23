#ifndef BLEND_H
#define BLEND_H

template<class T>
concept StatType = requires(T s, double v) {
    s.Update(v);
    { s.Get() };
};

template <StatType Stats>
class BlendExp
{
  static constexpr double EPS=1E-8;
  public:
    BlendExp(std::size_t n,double alpha,double beta)
    :n(n),beta(beta),px(0.0),
    x(n),w(n),zm(n),
    rsum(n,Stats(alpha))
    {
      if (n==0)
        throw std::invalid_argument("BlendExp: n must be >0");
      else
        std::fill(begin(w),end(w),1.0/n); //init equal weight
    };
    double Predict(const vec1D &input)
    {
      x=input;
      px=slmath::dot(x,w);
      return px;
    }
    void Update(double target)
    {
      UpdateScores(target);
      UpdateWeights();
      //UpdateBeta();
    }
    const vec1D &Weights()const {return w;}
  private:
    /*void UpdateBeta()
    {
      const double Hmax = std::log(n);
      const double Htarget = Hmax*0.3;
      double H = 0.0;
      for (std::size_t i=0;i<n;i++) {
        if (w[i]>1E-8) H-=w[i]*log(w[i]);
      }
      beta += 0.05*(H-Htarget);
      beta = std::clamp(beta,1.0,10.0);
    }*/
    void UpdateScores(double target)
    {
      //double loss_px = std::abs(target-px);

      for (std::size_t i=0;i<n;i++) {
        double loss_pi=std::abs(target-x[i]);

        // if score > 0 -> expert better then blend
        //double score=(loss_px - loss_pi);
        rsum[i].Update(-loss_pi);
      }
    }
    double calculate_z(const Stats &st) const
    {
      if constexpr (std::is_same_v<Stats, RunMeanVar>) {
        auto [mean,var] = st.Get();
        return beta*mean/(std::sqrt(var)+EPS);
      } else {
        return beta*st.Get();
      }
    }
    // softmax w_i = exp(-beta * normalized_regret)
    void UpdateWeights()
    {
      double max_z = -std::numeric_limits<double>::infinity();
      for (std::size_t i=0;i<n;i++) {
        zm[i] = calculate_z(rsum[i]);
        max_z = std::max(max_z,zm[i]);
      }

      //best expert has highest z-score -> weight=exp(0)=1
      double total=0.0;
      for (std::size_t i=0;i<n;i++) {
        w[i] = std::exp(zm[i]-max_z);
        total += w[i];
      }
      //normalize weights, total >= 1 from max-trick
      const double inv_total=1.0/total;
      for (double &val : w) {
          val *= inv_total;
      }
    }

    std::size_t n;
    double beta,px;
    vec1D x,w,zm;
    std::vector <Stats> rsum;
};

#endif

