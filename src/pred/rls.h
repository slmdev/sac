#ifndef RLS_H
#define RLS_H

#include "../global.h"
#include "../math/matrix.h"

// recursive least squares algorithm
class RLS {
  public:
    RLS(int n,double alpha);
    double Predict();
    double Predict(const std::vector<double> &pred);
    void Update(double val);
    void UpdateHist(double val);
  private:
    void UpdateGain();
    void UpdateP(); //update inverse of covariance matrix P(n)=1/lambda*P(n-1)-1/lambda * k(n)*x^T(n)*P(n-1)
    int n;
    double alpha,p;
    std::vector<double> hist,w,k;
    Matrix P;
};


#endif // RLS_H
