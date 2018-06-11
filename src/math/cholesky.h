#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "../global.h"
#include "../common/utils.h"
#include "matrix.h"

// calculate the cholesky decomposition of a symmetric pos. def. matrix
class Cholesky {
  public:
    Cholesky(int n)
    :n(n)
    {
      sol.resize(n);
      for (int i=0;i<n;i++) sol[i].resize(i+1);
      y.resize(n);
      regp=0.;
    };
    bool Test();
    int Factor(const Matrix &src);
    void Solve(const std::vector<double> &b,std::vector<double> &x);
    void Solve(const std::vector<double> &b);
    Matrix m;
    std::vector <std::vector<double>>sol;
    double regp;
  protected:
    void CalcY(const std::vector<double> &b);
    std::vector<double> y;
    double scale;
    int n;
};
#endif // CHOLESKY_H
