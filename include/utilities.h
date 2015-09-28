#ifndef UTILITIES_H_
#define UTILITIES_H_
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <valarray>
#include <assert.h>

static const char fINT[] = "%i";

template <class T>
int locate (const std::vector<T> &v, const T x){
  size_t n = v.size ();
  int jl = -1;
  int ju = n;
  bool as = (v[n-1] >= v[0]);
  while (ju-jl > 1){
    int jm = (ju+jl)/2;
    if ((x >= v[jm]) == as)
      jl=jm;
    else
      ju=jm;
  }
  if (x == v[0])
    return 0;
  else if (x == v[n-1])
    return n-2;
  else
    return jl;
}

template <class T> std:: string conv (T &val, const char *fact){
  char VAL[20]; sprintf (VAL, fact, val);
  return std:: string(VAL);
}

template <class T>
void fill_linear ( std::vector<T> &v, size_t n, T min, T max ){
  v.resize ( n );
  for ( size_t i = 0; i < n; i++ )
    v[i] = min + (max - min) * T (i)/T (n-1);
}

double getY(std:: vector<double> x, std:: vector<double> y,double xi);
float getYFloat(std:: vector<float> x, std:: vector<float> y,float xi);

#endif
