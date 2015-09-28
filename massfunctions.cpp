#include "massfunctions.h"
#include "utilities.h"

double st99nufnu(double nu, double a, double p, double a0){
  double nup = nu*a;
  //std:: cout << a << "   " << p << "  " << a0 << std:: endl;
  //exit(1);
  return a0*(1+1/pow(nup,p))*sqrt(nup/2.)*exp(-nup/2.)/sqrt(M_PI);
}
