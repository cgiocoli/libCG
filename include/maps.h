#ifndef MAPS_H_
#define MAPS_H_
#include <vector>
#include <valarray>
#include <string>
#include <stdio.h>     
#include "utilities.h"
#include <gsl/gsl_math.h>
#include <CCfits/CCfits>

using namespace CCfits;

void halve(std:: valarray<float> &imap, int inx, int iny, std:: valarray<float> &omap);

void writeFits(std:: string filename,std:: valarray<float> f, int npix, int npixy);

void readFits (std::string fn, std:: valarray<float> &map, int &nx, int &ny, double &fov);

void readFits (std::string fn, std:: valarray<float> &map, int &nx, int &ny, double &fov, double &cvir);

void readFits (std::string fn, std:: valarray<float> &map, int &nx, int &ny);

void stats(std:: valarray<float> map, float &mean, float &sigma, float &kurt, float &skew);

void stats_med(std:: vector<double> map, double &median, double &q25, double &q75,double &q025, double &q975);

void get_hist(std:: vector<float> map, int nbin, std:: vector<float> &xi, 
	      std:: vector<float> &yi);

#endif
