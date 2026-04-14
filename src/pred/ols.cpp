#include "ols.h"

//#define INIT_COV

//ordinary least squares using cholesky decomposition
//on a (weighted) covariance matrix estimate of input vectors
OLS::OLS(int n,int kmax,double lambda,double nu,double beta_sum,double beta_pow,double beta_add)
:x(n),
chol(n),
w(n),b(n),mcov(n,vec1D(n)),
n(n),kmax(kmax),lambda(lambda),nu((1.0-lambda)*nu),
beta_pow(beta_pow),beta_add(beta_add),esum(beta_sum),
pred(0.)
{
  km=0;
  #ifdef INIT_COV
    for (int i=0;i<n;i++) mcov[i][i]=1.0;
  #endif
}

double OLS::Predict()
{
  return (pred=slmath::dot(x,w));
}

void OLS::Update(double val)
{
  //running geometric sum of absolute prediction error
  const double e=val-pred;
  esum.Update(fabs(e));

  //generalized IRLS step, beta_pow=0->L2, beta_pow=1->L1
  const double c=pow(esum.Get()+beta_add,-beta_pow);
  const double ff=(1.0-lambda)*c;

  // update estimate of covariance matrix
  for (int j=0;j<n;j++) {
    const double xj=x[j];
    // only update lower triangular
    for (int i=0;i<=j;i++)
      mcov[j][i]=lambda*mcov[j][i]+ff*(xj*x[i]);

    b[j]=lambda*b[j]+ff*(xj*val);
  }
  km++;
  if (km>=kmax) {
    if (chol.Factor(mcov,nu))
      chol.Solve(b,w);
    km=0;
  }
}
