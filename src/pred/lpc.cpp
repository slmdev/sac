#include "lpc.h"

//#define INIT_COV

OLS::OLS(int n,int kmax,double lambda,double nu,double beta_sum,double beta_pow,double beta_add)
:x(n),
chol(n),
w(n),b(n),mcov(n,vec1D(n)),
n(n),kmax(kmax),lambda(lambda),nu(n*nu*(1.0-lambda)),
beta_pow(beta_pow),beta_add(beta_add),esum(beta_sum)
{
  km=0;
  pred=0.0;
  #ifdef INIT_COV
    for (int i=0;i<n;i++) mcov[i][i]=1.0;
  #endif
}

double OLS::Predict()
{
  pred=slmath::dot(x,w);
  return pred;
}

void OLS::Update(double val)
{
  //running geometric sum of absolute prediction error
  esum.Update(fabs(val-pred));
  //forgetting factor
  double ff=(1.0-lambda)*pow(esum.Get()+beta_add,-beta_pow);

  // update estimate of covariance matrix
  for (int j=0;j<n;j++) {
    // only update lower triangular
    for (int i=0;i<=j;i++) mcov[j][i]=lambda*mcov[j][i]+ff*(x[j]*x[i]);
    b[j]=lambda*b[j]+ff*(x[j]*val);
  }

  km++;
  if (km>=kmax) {
    if (!chol.Factor(mcov,nu)) chol.Solve(b,w);
    km=0;
  }
}
