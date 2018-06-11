#ifndef MATRIX_H
#define MATRIX_H

#include "../global.h"

// simple Square-Matrix class
class Matrix
{
  public:
    Matrix(){};
    Matrix(int dim):mat(dim,std::vector<double>(dim,0.))
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
    size_t GetDim()const {return mat.size();};
    std::vector<std::vector<double>> mat;
  private:
};

#endif // MATRIX_H
