#ifndef LPC_H
#define LPC_H

#include "../common/utils.h"

template <class T>
class OLS {
  typedef std::vector<T> vec1D;
  typedef std::vector<vec1D> vec2D;
  const T ftol=1E-8;
  public:
    OLS(int n,int kmax=1,T lambda=0.998,T nu=0.001,T beta_sum=0.6,T beta_pow=0.75,T beta_add=2)
    :x(n),n(n),kmax(kmax),lambda(lambda),nu(nu),
    w(n),b(n),mcov(n,vec1D(n)),mchol(n,vec1D(n)),
    beta_pow(beta_pow),beta_add(beta_add),esum(beta_sum)
    {
      km=0;
      pred=0.0;
    }
    /*T Predict(const vec1D &p) {
      x=p;
      pred=0.;
      for (int i=0;i<n;i++) pred+=w[i]*x[i];
      return pred;
    }*/
    T Predict() {
      pred=0.;
      for (int i=0;i<n;i++) pred+=w[i]*x[i];
      return pred;
    }
    void Update(T val)
    {
      T err=val-pred;
      esum.Update(fabs(err));
      double c0=pow(esum.Get()+beta_add,-beta_pow);

      for (int j=0;j<n;j++) {
        for (int i=0;i<n;i++) mcov[j][i]=lambda*mcov[j][i]+c0*(x[j]*x[i]);
        b[j]=lambda*b[j]+c0*(x[j]*val);
      }

      km++;
      if (km>=kmax) {
        if (!Factor(mcov)) Solve(b,w);
        km=0;
      }
    }
    vec1D x;
  //private:
    int Factor(const vec2D &mcov)
    {
      //std::cout << t << ' ' << nu << '\n';
      mchol=mcov; // copy the matrix
      for (int i=0;i<n;i++) mchol[i][i]+=nu;

      #if 1
      for (int i=0;i<n;i++) {
        for (int j=0;j<=i;j++) {
          T sum=mchol[i][j];
          for (int k=0;k<j;k++) sum-=(mchol[i][k]*mchol[j][k]);
          if (i==j) {
            if (sum>ftol) mchol[i][i]=sqrt(sum);
            else return 1;
          } else mchol[i][j]=sum/mchol[j][j];
        }
      }
      #else
      for (int i=0;i<n;i++) {
        for (int j=0;j<i;j++) {
          T sum=mchol[i][j];
          for (int k=0;k<j;k++) sum-=(mchol[i][k]*mchol[j][k]);
          mchol[i][j]=sum/mchol[j][j];
        }
        T sum=mchol[i][i];
        for (int k=0;k<i;k++) sum-=(mchol[i][k]*mchol[i][k]);
        if (sum>ftol) mchol[i][i]=sqrt(sum);
        else return 1;
      }
      #endif
      return 0;
    }

    void Solve(const vec1D &b,vec1D &x)
    {
      for (int i=0;i<n;i++) {
        T sum=b[i];
        for (int j=0;j<i;j++) sum-=(mchol[i][j]*x[j]);
        x[i]=sum/mchol[i][i];
      }
      for (int i=n-1;i>=0;i--) {
        T sum=x[i];
        for (int j=i+1;j<n;j++) sum-=(mchol[j][i]*x[j]);
        x[i]=sum/mchol[i][i];
      }
    }
    int n,kmax,km;
    T lambda,nu,pred;
    vec1D w,b;
    vec2D mcov,mchol;
    T beta_pow,beta_add;
    RunWeight esum;
};

#endif // LPC_H
