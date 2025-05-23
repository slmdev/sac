#ifndef MATH_H
#define MATH_H

#include "../global.h"
#include <cassert>

namespace slmath
{

  double dot_scalar(const vec1D &v1,const vec1D &v2)
  {
    assert(v1.size()==v2.size());
    double sum=0.0;
    for (std::size_t i=0;i<v1.size();++i)
      sum+=v1[i]*v2[i];
    return sum;
  }

  // vector = matrix * vector
  vec1D mul(const vec2D &m,const vec1D &v)
  {
    vec1D v_out(m.size());
    for (std::size_t i=0;i<m.size();i++)
      v_out[i]=slmath::dot_scalar(m[i],v);
    return v_out;
  }

  // vector = scalar * vector
  vec1D mul(const double s,const vec1D &v)
  {
    vec1D v_out(v.size());
    for (std::size_t i=0;i<v.size();i++)
      v_out[i]=s*v[i];
    return v_out;
  }

  // matrix = matrix  * matrix
  vec2D mul(const vec2D &m1, const vec2D &m2)
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

  // outer product of u*v^T
  vec2D outer(const vec1D &u,const vec1D &v)
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

