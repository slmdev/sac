#ifndef COV_H
#define COV_H

#include "matrix.h"

// adaptively calculate auto-covariance-matrix

// WLS is used to solved a weighted least squares problem
// covariance update R(n)=a*R(n-1)+x*x^T
// else the update of R is done via exp-smoothing i.e R(n)=a*R(n-1)+(1-a)x*x^T
#define WLS

class AutoCov {
  public:
    AutoCov(double alpha,int dim):c(dim),hist(dim),b(dim),alpha(alpha)
    {
      for (int i=0;i<dim;i++) c[i][i]=1.0;
    };
    void UpdateCov() {
      const int n=hist.size();
      for (int j=0;j<n;j++) {
        double *cj=&(c.mat[j][0]);
        for (int i=0;i<=j;i++) {
            #ifdef WLS
              cj[i]=alpha*cj[i]+(hist[i]*hist[j]);
            #else
              cj[i]=alpha*cj[i]+(1.0-alpha)*(hist[i]*histj);
            #endif
        }
      }
    }
    void UpdateB(double val)
    {
      const int n=hist.size();
      for (int i=0;i<n;i++) {
        #ifdef WLS
          b[i]=alpha*b[i]+(val*hist[i]);
        #else
          b[i]=alpha*b[i]+(1.0-alpha)*(val*hist[i]);
        #endif
          //b[i]+=0.001;
      }
    }
    void UpdateHist(double val,int n0,int n1) {
      for (int i=n1;i>n0;i--) hist[i]=hist[i-1];
      hist[n0]=val;
    }
    Matrix c;
    std::vector<double> hist,b;
    double alpha;
};
#endif // COV_H
