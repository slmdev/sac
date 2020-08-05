#ifndef GLOBAL_H
#define GLOBAL_H

//#include "windows.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <vector>

#define NDEBUG

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define EPS (1E-14)

#define DO(n) for (uint32_t _=0;_<n;_++)

typedef std::vector<double> vec1D;
typedef std::vector<std::vector<double>> vec2D;
typedef std::vector <double*>ptr_vec1D;

template <class T> T abs(T a){return ((a) < 0)?-(a):(a);}
template <class T> T div_signed(T val,T s){return val<0?-(((-val)+(1<<(s-1)))>>s):(val+(1<<(s-1)))>>s;};

#define RNDINT(f) ((int)(f >= 0.0 ? (f + 0.5))

#endif // GLOBAL_H
