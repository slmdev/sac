#ifndef RLS_H
#define RLS_H

#include "../global.h"

// simple Square-Matrix class
class SQMatrix
{
  public:
    SQMatrix(){};
    SQMatrix(int dim):mat(dim,std::vector<double>(dim,0.))
    {
    }
    std::vector<double>& operator[] (size_t i) { return mat[i]; } // index operator

    void Print() {
      int n=mat.size();
      for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) std::cout << mat[j][i] << " ";
        std::cout << std::endl;
      }
    }
    size_t Dim()const {return mat.size();};
    std::vector<std::vector<double>> mat;
  private:
};

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
    SQMatrix P;
};


#endif // RLS_H
