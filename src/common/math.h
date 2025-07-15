#ifndef MATH_H
#define MATH_H

#include "../global.h"
#include <cassert>
#include <cmath>
#include <immintrin.h>

namespace slmath
{

  // inplace cholesky
  // matrix must be positive definite and symmetric
  class Cholesky
  {
    public:
      const double ftol=1E-8;
      Cholesky(int n)
      :n(n),G(n,vec1D(n))
      {

      }
      int Factor(const vec2D &matrix,const double nu)
      {
        for (int i=0;i<n;i++) //copy lower triangular matrix
          std::copy_n(begin(matrix[i]),i+1,begin(G[i]));

        for (int i=0;i<n;i++) {

          // off-diagonal
          for (int j=0;j<i;j++) {
            double sum=G[i][j];
            for (int k=0;k<j;k++) sum-=(G[i][k]*G[j][k]);
            G[i][j]=sum/G[j][j];
          }

          // diagonal
          double sum=G[i][i]+nu; //add regularization
          for (int k=0;k<i;k++) sum-=(G[i][k]*G[i][k]);
          if (sum>ftol) G[i][i]=std::sqrt(sum);
          else return 1;
        }
        return 0;
      }
      void Solve(const vec1D &b,vec1D &x)
      {
        for (int i=0;i<n;i++) {
          double sum=b[i];
          for (int j=0;j<i;j++) sum-=(G[i][j]*x[j]);
          x[i]=sum/G[i][i];
        }
        for (int i=n-1;i>=0;i--) {
          double sum=x[i];
          for (int j=i+1;j<n;j++) sum-=(G[j][i]*x[j]);
          x[i]=sum/G[i][i];
        }
      }
      int n;
      vec2D G;
  };


  inline double dot(span_cf64 x,span_cf64 y)
  {
    assert(x.size()==y.size());
    const std::size_t n=x.size();
    double total=0.0;
    std::size_t i=0;

    if constexpr(SACGlobalCfg::USE_AVX2) {
      if constexpr(SACGlobalCfg::UNROLL_AVX2) {
        if (n>=8)
        {
          __m256d sum1 = _mm256_setzero_pd();
          __m256d sum2 = _mm256_setzero_pd();
          for (;i + 8 <= n;i += 8)
          {
            __m256d vx1 = _mm256_loadu_pd(&x[i]);
            __m256d vy1 = _mm256_loadu_pd(&y[i]);
            sum1 = _mm256_fmadd_pd(vx1, vy1, sum1);
            __m256d vx2 = _mm256_loadu_pd(&x[i + 4]);
            __m256d vy2 = _mm256_loadu_pd(&y[i + 4]);
            sum2 = _mm256_fmadd_pd(vx2, vy2, sum2);
          }
          sum1 = _mm256_add_pd(sum1, sum2);
          alignas(32) double buffer[4];
          _mm256_store_pd(buffer, sum1);
          total = buffer[0] + buffer[1] + buffer[2] + buffer[3];
        }
      } else if (n>=4)
      {
        __m256d sum = _mm256_setzero_pd();
        for (;i + 4 <= n;i += 4)
        {
          __m256d vx = _mm256_loadu_pd(&x[i]);
          __m256d vy = _mm256_loadu_pd(&y[i]);
          sum = _mm256_fmadd_pd(vx, vy, sum);
        }
        alignas(32) double buffer[4];
        _mm256_store_pd(buffer, sum);
        total = buffer[0] + buffer[1] + buffer[2] + buffer[3];
      }
    }

    for (;i<n;i++)
      total+=x[i]*y[i];
    return total;
  }

  // vector = matrix * vector
  inline vec1D mul(const vec2D &m,const vec1D &v)
  {
    vec1D v_out(m.size());
    for (std::size_t i=0;i<m.size();i++)
      v_out[i]=slmath::dot(m[i],v);
    return v_out;
  }

  // vector = scalar * vector
  inline vec1D mul(const double s,const vec1D &v)
  {
    vec1D v_out(v.size());
    for (std::size_t i=0;i<v.size();i++)
      v_out[i]=s*v[i];
    return v_out;
  }

  // matrix = matrix  * matrix
  inline vec2D mul(const vec2D &m1, const vec2D &m2)
  {
    vec2D m_out(m1.size(), vec1D(m2[0].size()));
    for (int j=0;j<(int)m_out.size();j++)
      for (int i=0;i<(int)m_out[0].size();i++)
      {
        double sum=0;
        for (int k=0;k<(int)m2.size();k++)
          sum += m1[j][k]*m2[k][i];
        m_out[j][i] = sum;
      }
    return m_out;
  }

  //vector s1*v1 + s2*v2
  inline vec1D mul_add(double s1,const vec1D &v1,double s2,const vec1D &v2)
  {
    assert(v1.size()==v2.size());
    vec1D v_out(v1.size());
    for (std::size_t i=0;i<v1.size();i++) {
      v_out[i] = s1*v1[i] + s2*v2[i];
    }
    return v_out;
  }

  // matrix s1*m1 + s2*m2
  inline vec2D mul_add(double s1,const vec2D &m1,double s2,const vec2D &m2)
  {
   assert(m1.size()==m2.size());
   vec2D m_out(m1.size());
   for (std::size_t j=0;j<m1.size();j++)
     m_out[j] = mul_add(s1,m1[j],s2,m2[j]);

   return m_out;
  }

  // outer product of u*v^T
  inline vec2D outer(const vec1D &u,const vec1D &v)
  {
    int nrows=u.size();
    int ncols=v.size();
    vec2D m_out(nrows, vec1D(ncols));
    for (int j=0;j<nrows;j++)
      for (int i=0;i<ncols;i++)
        m_out[j][i]=u[j]*v[i];
    return m_out;
  }

};

#endif

