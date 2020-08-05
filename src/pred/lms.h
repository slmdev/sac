#ifndef LMS_H
#define LMS_H

#include "../global.h"
#include "../common/histbuf.h"
#include <cmath>

class NLMS_ROLL {
  const double eps_pow=1.0;
  double sgn(double x) {
    if (x>0) return 1.;
    if (x<0) return -1.;
    return 0;
  }
  public:
    NLMS_ROLL(int n,double mu,double mu_decay=1.0,double pow_decay=0.8)
    :n(n),x(n),w(n),dw(n),mutab(n),powtab(n),mu(mu)
    {
      sum_powtab=0;
      for (int i=0;i<n;i++) {
         powtab[i]=1.0/(pow(1+i,pow_decay));
         sum_powtab+=powtab[i];
         mutab[i]=pow(mu_decay,i);
      }
    }
    double Predict()
    {
      pred=0.;
      for (int i=0;i<n;i++) pred+=w[i]*x[i];
      return pred;
    }
    void Update(double val) {
      double spow=0.0;
      for (int i=0;i<n;i++) {
        spow+=powtab[i]*(x[i]*x[i]);
      }
      const double wgrad=(val-pred)*sum_powtab/(eps_pow+spow);
      for (int i=0;i<n;i++) {
        dw[i]=(mu*mutab[i]*wgrad*x[i]); //+mom*dw[i];
        w[i]+=dw[i];
      }
      x.PushBack(val);
    };
    ~NLMS_ROLL(){};
  protected:
    int n;
    RollBuffer <double>x;
    vec1D w,dw,mutab,powtab;
    double pred,sum_powtab;
    double mu;
};

class LMS_ROLL {
  double sgn(double x) {
    if (x>0) return 1.;
    else if (x<0) return -1;
    else return 0;
  }
  const double eps_pow=1.0;
  const double eps_ada=1.0;
  public:
    LMS_ROLL(int n,double mu,double beta=0.95,double mu_decay=1.0,double pow_decay=0.8)
    :n(n),x(n),w(n),eg(n),mutab(n),powtab(n),mu(mu),beta(beta)
    {
      do_lms=false;
      do_ada=false;

      sum_powtab=0;
      for (int i=0;i<n;i++) {
         powtab[i]=1.0/(pow(1+i,pow_decay));
         sum_powtab+=powtab[i];
         mutab[i]=pow(mu_decay,i);
      }
    }
    double Predict()
    {
      pred=0.0;
      for (int i=0;i<n;i++) pred+=w[i]*x[i];
      return pred;
    }
    void Update(double val) {
      const double err=val-pred; // prediction error

      #if 1
      double wgrad=0.0;
      if (do_lms) { // lms
        wgrad=err;
      } else { // nlms
        double spow=0.0;
        for (int i=0;i<n;i++) {
          spow+=powtab[i]*(x[i]*x[i]);
        }
        wgrad=err*sum_powtab/(n*eps_pow+spow);
      }

      if (do_ada) {
        for (int i=0;i<n;i++) {
          //double const grad=sgn(wgrad)*x[i];//-0.001*sgn(w[i]);
          double const grad=wgrad*x[i];//-0.001*sgn(w[i]);
          eg[i]=beta*eg[i]+(1.0-beta)*grad*grad; //accumulate squared gradients
          w[i]+=mu*mutab[i]*grad/(sqrt(eg[i])+eps_ada);
        }
      } else {
        for (int i=0;i<n;i++) {
          w[i]+=mu*mutab[i]*wgrad*x[i];
        }
      }
      #else
      for (int i=0;i<n;i++) {
        double const grad=err*x[i]; // gradient + l1-regularization
        eg[i]=beta*eg[i]+(1.0-beta)*grad*grad; //accumulate squared gradients
        double g=grad*1.0/(sqrt(eg[i])+eps_ada);// update weights
        w[i]+=mu*g;
      }
      #endif
      x.PushBack(val);
    };
    ~LMS_ROLL(){};
  protected:
    bool do_lms,do_ada;
    int n;
    RollBuffer <double>x;
    vec1D w,eg,mutab,powtab;
    double pred,sum_powtab;
    double mu,beta;
};

class LMS {
  protected:
    inline double sgn(double x) {
      if (x>0) return 1.;
      if (x<0) return -1;
      return 0;
    }
  public:
    LMS(int n,double mu)
    :n(n),x(n),w(n),mu(mu)
    {
    }
    virtual double Predict(const vec1D &inp)
    {
      x=inp;
      pred=0.0;
      for (int i=0;i<n;i++) pred+=w[i]*x[i];
      return pred;
    }
    virtual void Update(double)=0;
    virtual ~LMS(){};
    int n;
    vec1D x,w;
  protected:
    double mu,pred;
};

class LMS_ADA : public LMS
{
  public:
    LMS_ADA(int n,double mu,double beta=0.95,double nu=0.0)
    :LMS(n,mu),eg(n),beta(beta),nu(nu)
    {
    }
    virtual void Update(double val) {
      const double err=val-pred; // prediction error
      for (int i=0;i<n;i++) {
        double const grad=err*x[i]-nu*sgn(w[i]); // gradient + l1-regularization
        eg[i]=beta*eg[i]+(1.0-beta)*grad*grad; //accumulate gradients
        double g=grad*1.0/(sqrt(eg[i])+1E-5);// update weights
        w[i]+=mu*g;
      }
    }
  protected:
    vec1D eg;
    double beta,nu;
};

class LAD_ADA : public LMS
{
  public:
    LAD_ADA(int n,double mu,double beta=0.95)
    :LMS(n,mu),eg(n),beta(beta)
    {
    }
    virtual void Update(double val) {
      const double err=val-pred; // prediction error
      for (int i=0;i<n;i++) {
        double const grad=sgn(err)*x[i];
        eg[i]=beta*eg[i]+(1.0-beta)*grad*grad; //accumulate gradients
        double g=grad*1.0/(sqrt(eg[i])+1E-5);// update weights
        w[i]+=mu*g;
      }
    }
  protected:
    vec1D eg;
    double beta;
};

class LMS_ADAM : public LMS
{
  public:
    LMS_ADAM(int n,double mu,double beta1=0.9,double beta2=0.999)
    :LMS(n,mu),M(n),S(n),beta1(beta1),beta2(beta2)
    {
      power_beta1=1.0;
      power_beta11=beta1;
      power_beta2=1.0;
    }
    void Update(double val) {
      power_beta1*=beta1;
      power_beta11*=beta1;
      power_beta2*=beta2;
      const double err=val-pred; // prediction error
      for (int i=0;i<n;i++) {
        double const grad=err*x[i]; // gradient

        M[i]=beta1*M[i]+(1.0-beta1)*grad;
        S[i]=beta2*S[i]+(1.0-beta2)*(grad*grad);

        /*double m_hat=beta1*M[i]/(1.0-power_beta11)+((1.0-beta1)*grad/(1.0-power_beta1));
        double n_hat=beta2*S[i]/(1.0-power_beta2);*/
        double m_hat=M[i]/(1.0-power_beta1);
        double n_hat=S[i]/(1.0-power_beta2);
        w[i]+=mu*m_hat/(sqrt(n_hat)+1E-5);
      }
    }
  private:
    vec1D M,S;
    double beta1,beta2,power_beta1,power_beta11,power_beta2;
};

// sign-sign lms algorithm
class SSLMS : public LMS {
  public:
      SSLMS(int n,double mu)
      :LMS(n,mu)
      {
      }
      void Update(double val)
      {
        double e=val-pred;
        const double wf=mu*sgn(e);
        for (int i=0;i<n;i++) {
           w[i]+=wf*sgn(x[i]);
        }
      }
};

#endif // LMS_H
