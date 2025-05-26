#include "rls.h"
#include "../common/math.h"
#include "../common/utils.h"

RLS::RLS(int n,double gamma,double nu)
:n(n),
px(0.),
hist(n),w(n),
P(n,vec1D(n)), // inverse covariance matrix
alc(gamma)
{
  for (int i=0;i<n;i++)
    P[i][i]=1.0/nu;
}

double RLS::Predict()
{
  px=slmath::dot_scalar(hist,w);
  return px;
}

double RLS::Predict(const vec1D &input)
{
  hist=input;
  return Predict();
}


void RLS::Update(double val)
{
  double err=val-px;

  vec1D vt1=slmath::mul(P,hist); //phi=hist P hist
  // a priori variance of prediction
  double phi=slmath::dot_scalar(hist,vt1);

  // Normalized Innovation Squared
  // quantifies how "unexpected" the observation is
  // relative to the current uncertainty
  double metric = (err*err);///(phi+1E-5);
  double alpha=alc.update(metric);

  //update inverse of covariance matrix
  //P(n)=1/lambda*P(n-1)-1/lambda * k(n)*x^T(n)*P(n-1)
  double denom=1./(alpha+phi);
  double inv_alpha=1.0/(alpha);
  for (int i=0;i<n;i++)
    for (int j=0;j<=i;j++) {
      double m=vt1[i]*vt1[j]; // outer product of vt1
      double v=(P[i][j] - denom * m) * inv_alpha;
      P[i][j] = v;
      P[j][i] = v;
    }

  // update weights
  for (int i=0;i<n;i++)
      w[i]+=err*(denom*vt1[i]);
}

void RLS::UpdateHist(double val)
{
  Update(val);
  miscUtils::RollBack(hist,val);
}
