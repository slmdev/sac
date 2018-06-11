#include "cholesky.h"

bool Cholesky::Test()
{
  Matrix mtmp;
  mtmp.mat={{9,3,-6,12},{3,26,-7,-11},{-6,-7,9,7},{12,-11,7,65}};
  std::vector<double> b={72,34,22,326};
  std::vector<double> sol={2,4,3,5};
  std::vector<double> x(4);
  Factor(mtmp);
  Solve(b,x);
  if (MathUtils::SumDiff(sol,x)<EPS) return true;
  else return false;
}

int Cholesky::Factor(const Matrix &msrc)
{
  int i,j,k;
  m=msrc; // copy the matrix
  for (int i=0;i<n;i++) m[i][i]+=regp;

  double sum;
  for (i=0;i<n;i++) {
    for (j=0;j<i;j++) {
      sum=m[i][j];
      for (k=0;k<j;k++) sum-=(m[i][k]*m[j][k]);
      m[i][j]=sum/m[j][j];
    }
    sum=m[i][i];
    for (k=0;k<i;k++) sum-=(m[i][k]*m[i][k]);
    if (sum>0.) m[i][i]=sqrt(sum);
    else return 1; // matrix indefinit
  }
  return 0;
}

void Cholesky::CalcY(const std::vector<double> &b)
{
  int i,j;
  double sum;
  for (i=0;i<n;i++) {
     sum=b[i];
     for (j=0;j<i;j++) sum-=(m[i][j]*y[j]);
     y[i]=sum/m[i][i];
  }
}

void Cholesky::Solve(const std::vector<double> &b,std::vector<double> &x)
{
  //vector<double> bx=b;
  //for (int i=0;i<b.size();i++) bx[i]=bx[i]+regp/double(n);
  CalcY(b);
  int i,j;
  double sum;
  for (i=n-1;i>=0;i--) {
    sum=y[i];
    for (j=i+1;j<n;j++) sum-=(m[j][i]*x[j]);
    x[i]=sum/m[i][i];
  }
}

// solve for all orders
void Cholesky::Solve(const std::vector<double> &b)
{
  CalcY(b);

  int i,j;
  double sum;
  for (int o=0;o<n;o++)
  {
    for (i=o;i>=0;i--) {
      sum=y[i];
      for (j=i+1;j<=o;j++) sum-=(m[j][i]*sol[o][j]);
      sol[o][i]=sum/m[i][i];
    }
  }
}
