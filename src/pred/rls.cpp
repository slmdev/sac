#include "rls.h"

double VecMulVecInner(const std::vector<double> &v1,const std::vector<double> &v2)
{
  double sum=0.;
  for (size_t i=0;i<v1.size();i++) sum+=v1[i]*v2[i];
  return sum;
}

std::vector<double> MatrixMulVec(Matrix &m,const std::vector<double> &v)
{
  int n=v.size();
  std::vector<double> t(n);
  for (int i=0;i<n;i++) t[i]=VecMulVecInner(m[i],v);
  return t;
}

std::vector<double> VecMulMatrix(const std::vector<double> &v, Matrix &m)
{
  int n=v.size();
  std::vector<double> t(n);
  for (int j=0;j<n;j++) {
    double sum=0;
    for (int i=0;i<n;i++) sum+=v[i]*m[i][j];
    t[j]=sum;
  }
  return t;
}

std::vector<double> ScalarMulVec(const double x,const std::vector<double> &v)
{
  int n=v.size();
  std::vector<double> t(n);
  for (int i=0;i<n;i++) t[i]=x*v[i];
  return t;
}

Matrix VecMulVecOuter(const std::vector<double> &v1,const std::vector<double> &v2)
{
  int n=v1.size();
  Matrix m(n);
  for (int j=0;j<n;j++)
    for (int i=0;i<n;i++) m[j][i]=v1[j]*v2[i];
  return m;
}

Matrix MatrixMulMatrix(Matrix &m1,Matrix &m2)
{
  int n=m1[0].size();
  Matrix m(n);
  for (int i=0;i<n;i++) { // for every row
    for (int j=0;j<n;j++) { // for every column
      double sum=0.;
      for (int k=0;k<n;k++) sum+=m1[i][k]*m2[k][j]; // multiply row*column
      m[i][j]=sum;
    }
  }
  return m;
}

RLS::RLS(int n,double alpha)
:n(n),alpha(alpha),hist(n),w(n),P(n)
{
  for (int i=0;i<n;i++) P[i][i]=1;
}

double RLS::Predict()
{
  p=0.;
  for (int i=0;i<n;i++) p+=hist[i]*w[i];
  return p;
}

double RLS::Predict(const std::vector<double> &pred)
{
  hist=pred;
  return Predict();
}

void RLS::UpdateGain()
{
 #if 0
  vector<double> vt1=MatrixMulVec(P,hist); // P(n-1)*x(n)
  vector<double> vt2=VecMulMatrix(hist,P); // x^T(n)*P(n-1)
  double t=VecMulVecInner(vt2,hist); // x^T(n)*P(n-1)*x(n)
  k=ScalarMulVec(1./(alpha+t),vt1);
 #else
  std::vector<double> vt1=MatrixMulVec(P,hist);
  k=ScalarMulVec(1./(alpha+VecMulVecInner(hist,vt1)),vt1);
 #endif
}

      //update inverse of covariance matrix P(n)=1/lambda*P(n-1)-1/lambda * k(n)*x^T(n)*P(n-1)
void RLS::UpdateP()
{
 Matrix m1=VecMulVecOuter(k,hist); //m1 is symmetric
 Matrix m2=MatrixMulMatrix(m1,P);
 for (int i=0;i<n;i++)
 for (int j=0;j<n;j++) P[i][j]=1./alpha*(P[i][j]-m2[i][j]);
}

void RLS::Update(double val)
{
  UpdateGain();
  UpdateP();

  for (int i=0;i<n;i++) w[i]=w[i]+(val-p)*k[i]; //apply gain to the weights
}

void RLS::UpdateHist(double val)
{
  Update(val);
  for (int i=n-1;i>0;i--) hist[i]=hist[i-1];hist[0]=val;
}
