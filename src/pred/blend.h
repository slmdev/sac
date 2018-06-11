#ifndef BLEND_H
#define BLEND_H

#include "lpc.h"


class BlendErr {
  public:
      BlendErr(int n,double alpha)
      :serr(n+1,RunExp(alpha)),w(n),p(n),n(n)
      {
      }
      double Predict(std::vector<double> pred)
      {
        p=pred;
        pm=0.;
        for (int i=0;i<n;i++) pm+=w[i]*p[i];

        int mini=0;
        double minc=serr[0].Get();
        for (int i=1;i<=n;i++) {
            double t=serr[i].Get();
            if (t<minc) {
                minc=t;
                mini=i;
            }
        }
        return pm;
      }
      void Update(int32_t val)
      {
        for (int i=0;i<n;i++) {
          double e=(val-p[i]);
          serr[i].Update(e*e);
        }
        serr[n].Update((val-pm)*(val-pm));

        double wsum=0;
        for (int i=0;i<n;i++) {
          double e=serr[i].Get();
          double t=10000./(e+0.0001);
          w[i]=t*t;
          wsum+=w[i];
        }
        if (wsum>0) for (int i=0;i<n;i++) w[i]/=wsum;
      }
  protected:
    std::vector <RunExp> serr;
    std::vector<double> w,p;
    double pm;
    int n;

};

class BlendRPROP {
    inline double sgn(double v){
      if (v<0.) return -1;
      if (v>0.) return 1;
      return 0;
    }
  public:
      BlendRPROP(int n,double mu,double mu_min,double mu_max)
      :w(n,0.0),y(n,mu),p(n),lgrad(n),mu(mu),mu_min(mu_min),mu_max(mu_max),n(n)
      {
      }
      double Predict(const std::vector<double>&pred)
      {
        p=pred;
        pm=0.;
        for (int i=0;i<n;i++) pm+=w[i]*p[i];
        return pm;
      }
      void Update(double val)
      {
        double e=val-pm;
        for (int i=0;i<n;i++) {
            const double grad=e*p[i];
            int sg=grad*lgrad[i];
            if (sg>0) {
              y[i]=std::min(y[i]*1.2,mu_max);
            } else if (sg<0) {
              y[i]=std::max(y[i]*0.5,mu_min);
            }
            w[i]+=y[i]*sgn(grad);
            lgrad[i]=grad;
        }
      }
  private:
    std::vector<double> w,y;
    std::vector<double> p,lgrad;
    double mu,mu_min,mu_max,pm;
    int32_t n;
};

class BlendRMS {
  public:
      BlendRMS(int n,double mu)
      :w(n,0.0),eg(n),p(n),mu(mu),n(n)
      {

      }
      int32_t Predict(const std::vector<int32_t>&pred)
      {
        p=pred;
        double sum=0.;
        for (int i=0;i<n;i++) sum+=w[i]*p[i];
        pm=round(sum);
        return pm;
      }
      void Update(int32_t val)
      {
        const double err=val-pm;
        const double rho=0.9;
        for (int i=0;i<n;i++) {
          double const grad=err*p[i]; // gradient
          eg[i]=rho*eg[i]+(1.0-rho)*grad*grad; //accumulate gradients
          w[i]+=(mu*grad)/sqrt(eg[i]+0.1);// update weights
        }
      }
  private:
    std::vector<double> w,eg;
    std::vector <int32_t>p;
    double mu;
    int32_t n,pm;
};


// sign-sign lms algorithm
class SSLMS {
    inline double sgn(double v){
      if (v<0) return -1;
      if (v>0) return 1;
      return 0;
    }
  public:
      SSLMS(int n,double mu)
      :w(n,0.0),n(n),mu(mu),p(n)
      {
      }
      double Predict(const std::vector<double> &pred) {
        p=pred;
        pm=0.;
        for (int i=0;i<n;i++) pm+=w[i]*p[i];
        return pm;
      }
      void Update(double val)
      {
        double e=val-pm;
        const double wf=mu*sgn(e);
        for (int i=0;i<n;i++) {
           w[i]+=wf*sgn(p[i]);
        }
      }
    std::vector<double> w;
  private:
    int n;
    double mu,pm;
    std::vector<double> p;
    //RunWeight rwa,rwb;
};


#endif // BLEND_H
