#include <cmath>
#include <iostream>
#include <algorithm>
#include "utilities.h"

double getY(std:: vector<double> x, std:: vector<double> y,double xi){
  int nn = x.size();
  if(x[0]<x[nn-1]){
    if(xi>x[nn-1]) return y[nn-1];
    if(xi<x[0]) return y[0];
  }
  else{
    if(xi<x[nn-1]) return y[nn-1];
    if(xi>x[0]) return y[0];
  }
  int i = locate (x,xi);
  i = std::min (std::max (i,0), int (nn)-2);
  double f=(xi-x[i])/(x[i+1]-x[i]);
  if(i>1 && i<nn-2){
    double a0,a1,a2,a3,f2;
    f2 = f*f;
    a0 = y[i+2] - y[i+1] - y[i-1] + y[i];
    a1 = y[i-1] - y[i] - a0;
    a2 = y[i+1] - y[i-1];
    a3 = y[i];
    return a0*f*f2+a1*f2+a2*f+a3;
  }
  else return f*y[i+1]+(1-f)*y[i];
}

float getYFloat(std:: vector<float> x, std:: vector<float> y,float xi){
  int nn = x.size();
  if(x[0]<x[nn-1]){
    if(xi>x[nn-1]) return y[nn-1];
    if(xi<x[0]) return y[0];
  }
  else{
    if(xi<x[nn-1]) return y[nn-1];
    if(xi>x[0]) return y[0];
  }
  int i = locate (x,xi);
  i = std::min (std::max (i,0), int (nn)-2);
  float f=(xi-x[i])/(x[i+1]-x[i]);
  if(i>1 && i<nn-2){
    float a0,a1,a2,a3,f2;
    f2 = f*f;
    a0 = y[i+2] - y[i+1] - y[i-1] + y[i];
    a1 = y[i-1] - y[i] - a0;
    a2 = y[i+1] - y[i-1];
    a3 = y[i];
    return a0*f*f2+a1*f2+a2*f+a3;
  }
  else return f*y[i+1]+(1-f)*y[i];
}
