#ifndef LMS_H
#define LMS_H

#include <cmath>
#include "../global.h"
#include "../common/histbuf.h"
#include "../common/utils.h"

class LS_Stream {
  public:
    LS_Stream(int n)
    :n(n),x(n),w(n),pred(0.)
    {

    }
    double Predict()
    {
      pred=MathUtils::dot(x.data(),w.data(),n);
      return pred;
    }
    virtual void Update(double val)=0;
    virtual ~LS_Stream(){};
  protected:
    int n;
    RollBuffer2<double>x;
    std::vector<double,align_alloc<double>> w;
    double pred;
};

/*
powtab[0] = 1.0;
double r = pow(n, -pow_decay / double(n-1));
for (int i = 1; i < n; i++) {
    powtab[i] = powtab[i-1] * r;
}

void update_w_avx(double* w, const double* mutab, const double* x, double wgrad, int n)
{
    __m256d wgrad_vec = _mm256_set1_pd(wgrad);
    int i = 0;


    for (; i <= n - 4; i += 4) {
        __m256d w_vec = _mm256_loadu_pd(&w[i]);
        __m256d mutab_vec = _mm256_loadu_pd(&mutab[i]);
        __m256d x_vec = _mm256_loadu_pd(&x[i]);

        __m256d wx = _mm256_mul_pd(wgrad_vec, x_vec);
        __m256d result = _mm256_fmadd_pd(mutab_vec, wx, w_vec);
        _mm256_storeu_pd(&w[i], result);
    }

    for (; i < n; ++i) {
        w[i] += mutab[i] * (wgrad * x[i]);
    }
}*/


class NLMS_Stream : public LS_Stream
{
  const double eps_pow=1.0;
  public:
    NLMS_Stream(int n,double mu,double mu_decay=1.0,double pow_decay=0.8)
    :LS_Stream(n),mutab(n),powtab(n),mu(mu)
    {
      sum_powtab=0;
      for (int i=0;i<n;i++) {
         powtab[i]=1.0/(pow(1+i,pow_decay));
         sum_powtab+=powtab[i];
         mutab[i]=pow(mu_decay,i);
      }
    }

  #if defined(USE_AVX256)
    double calc_spow(const double *x,const double *powtab,std::size_t n)
    {
      double spow=0.0;
      std::size_t i=0;
      if (n>=4) {
        __m256d sum_vec = _mm256_setzero_pd();
        for (; i + 4 <= n; i += 4) {
          __m256d x_vec = _mm256_loadu_pd(&x[i]);
          __m256d pow_vec = _mm256_load_pd(&powtab[i]);
          __m256d x_squared = _mm256_mul_pd(x_vec, x_vec);
          sum_vec = _mm256_fmadd_pd(pow_vec, x_squared, sum_vec);
        }

        alignas(32) double buffer[4];
        _mm256_store_pd(buffer, sum_vec);
        spow = buffer[0] + buffer[1] + buffer[2] + buffer[3];
      }

      for (;i<n;i++)
        spow += powtab[i] * (x[i] * x[i]);
      return spow;
    }
  #else
    double calc_spow(const double *x,const double *powtab,std::size_t n)
    {
      double spow=0.0;
      for (std::size_t i=0;i<n;i++) {
        spow+=powtab[i]*(x[i]*x[i]);
      }
      return spow;
    }
  #endif

    void Update(double val) override
    {
      const double spow=calc_spow(x.data(),powtab.data(),n);
      const double wgrad=mu*(val-pred)*sum_powtab/(eps_pow+spow);
      for (int i=0;i<n;i++) {
        w[i]+=mutab[i]*(wgrad*x[i]);
      }
      x.push(val);
    };
    ~NLMS_Stream(){};
  protected:
    std::vector<double,align_alloc<double>> mutab,powtab;
    double sum_powtab;
    double mu;
};

class LADADA_Stream : public LS_Stream
{
  public:
    LADADA_Stream(int n,double mu,double beta=0.97)
    :LS_Stream(n),eg(n),mu(mu),beta(beta)
    {

    }
    void Update(double val) override
    {
      const double serr=MathUtils::sgn(val-pred); // prediction error
      for (int i=0;i<n;i++) {
        double const grad=serr*x[i];
        eg[i]=beta*eg[i]+(1.0-beta)*grad*grad; //accumulate gradients
        double g=grad*1.0/(sqrt(eg[i])+1E-5);// update weights
        w[i]+=mu*g;
      }
      x.push(val);
    }
  protected:
    vec1D eg;
    double mu,beta;
};

class LMSADA_Stream : public LS_Stream
{
  public:
    LMSADA_Stream(int n,double mu,double beta=0.97,double nu=0.0)
    :LS_Stream(n),eg(n),mu(mu),beta(beta),nu(nu)
    {

    }
    void Update(double val) override
    {
      const double err=val-pred; // prediction error
      for (int i=0;i<n;i++) {
        double const grad=err*x[i]-nu*MathUtils::sgn(w[i]);
        eg[i]=beta*eg[i]+(1.0-beta)*grad*grad; //accumulate gradients
        double g=grad*1.0/(sqrt(eg[i])+1E-5);// update weights
        w[i]+=mu*g;
      }
      x.push(val);
    }
  protected:
    vec1D eg;
    double mu,beta,nu;
};


class LMS {
  protected:

  public:
    LMS(int n,double mu)
    :n(n),x(n),w(n),mu(mu),pred(0)
    {
    }
    double Predict(const vec1D &inp)
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
        double const grad=err*x[i]-nu*MathUtils::sgn(w[i]); // gradient + l1-regularization
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
      const double serr=MathUtils::sgn(val-pred); // prediction error
      for (int i=0;i<n;i++) {
        double const grad=serr*x[i];
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
        const double wf=mu*MathUtils::sgn(e);
        for (int i=0;i<n;i++) {
           w[i]+=wf*MathUtils::sgn(x[i]);
        }
      }
};

#endif // LMS_H
