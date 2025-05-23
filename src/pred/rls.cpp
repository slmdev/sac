#include "rls.h"
#include "../common/math.h"
#include "../common/utils.h"

 //update the sensitivity state zeta
  /*double denom = alpha + phi;
  for(int i=0;i<n;++i){
    double dk = -k[i] / denom;
    zeta[i] += err * dk;
  }

  // 4) update alpha via LMS on the meta gradient
  double hTz = slmath::dot_scalar(hist, zeta);
  alpha += alpha_mu * err * hTz;
  alpha = std::clamp(alpha,0.99,0.999);*/
  //std::cout << alpha << ' ';


  /*double S = phi + 1.0/alpha;
  double grad=(err*err - S) / (2.0 * alpha*alpha * S*S);
  alpha = alpha - 0.005*grad;
  alpha = std::clamp(alpha,0.99,0.999);
  //std::cout << alpha << ' ';*/
/*void RLS::UpdateP(double alpha,const vec1D &k)
{
 vec2D m1=slmath::outer(k,hist); //m1 is symmetric
 vec2D m2=slmath::mul(m1,P);

 for (int i=0;i<n;i++)
  for (int j=0;j<n;j++)
    P[i][j]=1./alpha*(P[i][j]-m2[i][j]);
}*/


RLS::RLS(int n,double alpha_mu,double nu)
:n(n),
px(0.),alpha_mu(alpha_mu),
hist(n),w(n),
P(n,vec1D(n)), // inverse covariance matrix
alc(alpha_mu)
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
  double metric = (err*err);///(phi+1.0);
  double alpha=alc.update(metric);

  //update inverse of covariance matrix
  //P(n)=1/lambda*P(n-1)-1/lambda * k(n)*x^T(n)*P(n-1)
  double denom=1./(alpha+phi);
  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++) {
      double m=vt1[i]*vt1[j]; // outer product of vt1
      P[i][j] = (P[i][j] - denom * m) / alpha;
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
