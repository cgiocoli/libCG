#include "maps.h"

void halve(std:: valarray<float> &imap, int inx, int iny, std:: valarray<float> &omap){
  if(inx % 2){
    // it is odd
    std:: cout << " nx of the input map is odd " << std:: endl;
    std:: cout << " I cannot process ... I will STOP here!!! " << std:: endl;
    exit(1);
  }
  if(iny % 2){
    // it is odd
    std:: cout << " ny of the input map is odd " << std:: endl;
    std:: cout << " I cannot process ... I will STOP here!!! " << std:: endl;
    exit(1);
  }

  int nx = inx/2;
  int ny = iny/2;

  omap.resize(nx*ny);
  for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){
      omap[i+nx*j] = (imap[(2*i)+inx*(2*j)] + 
		      imap[(2*i+1)+inx*(2*j)] + 
		      imap[(2*i)+inx*(2*j+1)] + 
		      imap[(2*i+1)+inx*(2*j+1)])/4;
      
    }
}

void writeFits(std:: string filename,std:: valarray<float> f, int npix, int npixy){
  long naxis=2;
  long naxes[2]={npix,npixy};
  std::auto_ptr<FITS> fout(new FITS(filename,FLOAT_IMG,naxis,naxes));
  std::vector<long> naxex(2);
  naxex[0]=npix;
  naxex[1]=npixy;
  PHDU *phout=&fout->pHDU();
  phout->write( 1, npix*npixy, f );
}

void readFits (std::string fn, std:: valarray<float> &map, int &nx, int &ny){
  std::auto_ptr<FITS> ff(new FITS (fn, Read));
  PHDU *h0=&ff->pHDU();
  nx=h0->axis(0);
  ny=h0->axis(1);
  h0->read(map);
}

// try to read also the field of view
void readFits (std::string fn, std:: valarray<float> &map, int &nx, int &ny, double &fov){
  std::auto_ptr<FITS> ff(new FITS (fn, Read));
  PHDU *h0=&ff->pHDU();
  nx=h0->axis(0);
  ny=h0->axis(1);
  h0->read(map);
  try {
    h0->readKey ("SIDEL",fov);
  }
  catch(CCfits::HDU::NoSuchKeyword) {
    fov=-1;
  }
}

// try to read also the field of view and the concentration
void readFits (std::string fn, std:: valarray<float> &map, int &nx, int &ny, double &fov, double &cvir){
  std::auto_ptr<FITS> ff(new FITS (fn, Read));
  PHDU *h0=&ff->pHDU();
  nx=h0->axis(0);
  ny=h0->axis(1);
  h0->read(map);
  try {
    h0->readKey ("SIDEL",fov);
  }
  catch(CCfits::HDU::NoSuchKeyword) {
    fov=-1;
  }
  try {
    h0->readKey ("HIERARCH CONCENTRATION",cvir);
  }
  catch(CCfits::HDU::NoSuchKeyword) {
    cvir=-1;
  }
}

// ... vector has to be sorted!
void stats_med(std:: vector<double> map, double &median, double &q25, double &q75,double &q025, double &q975){
  int n = map.size();
  if(n % 2){ 
    /* n is odd */ 
    int m = (n-1)/2;
    median = map[m];
    int m1 = (n-1)/4;
    q25 = map[m1];
    int m2 = 3*(n-1)/4;
    q75 = map[m2];
  }else{
    int m = n/2;
    median = 0.5*(map[m-1]+map[m]);
    int m1 = n/4;
    q25 = 0.5*(map[m1-1]+map[m1]);
    int m2 = 3*n/4;
    q75 = 0.5*(map[m2-1]+map[m2]);
  }
  int m025 = int(2.5/100.*(n-1));
  q025 = map[m025];
  int m975 = int(97.5/100.*(n-1));
  q975 = map[m975];
}

void stats(std:: valarray<float> map, float &mean, float &sigma, float &kurt, float &skew){
  float sum = map.sum();
  int n = map.size();
  mean = map.sum()/float(n);
  std:: valarray <float> maps(n);
  valarray<float> maps2(n),maps3(n),maps4(n);
  for(int i=0;i<n;i++){
    maps2[i] = gsl_pow_2(map[i] - mean);
    maps3[i] = gsl_pow_3(map[i] - mean);
    maps4[i] = gsl_pow_4(map[i] - mean);
  }
  sum = maps2.sum();
  sigma = sqrt(sum/(float(n)-1.));
  sum = maps3.sum();
  double mu3 = sum/(float(n)-1.);
  sum = maps4.sum();
  double mu4 = sum/(float(n)-1.);
  kurt = mu4/gsl_pow_4(sigma) -3;
  skew = mu3/gsl_pow_3(sigma);
}

void get_hist(std:: vector<float> map, int nbin, std:: vector<float> &xi, 
	      std:: vector<float> &yi){
  int n = map.size();
  std:: vector<float> extxi;
  // select min and max val !!!!!!!!!!!!!!!
  float min = map[0];
  float max = map[n-1];
  if(min==-1e+6){
    min = -1;
    max = 1.5;
  }
  fill_linear(extxi,nbin,min,max);
  float dx = (max-min)/(nbin-1);
  xi.resize(nbin-1);
  yi.resize(nbin-1);
  for(int l=0;l<nbin-1;l++){
    xi[l] = 0;
    yi[l] = 0;
    for(int i=0;i<n;i++){
      if(map[i]>=extxi[l] && map[i]<extxi[l+1]){
	yi[l]++;
	xi[l]+=map[i];
      }
      // since they are sorted
      if(map[i]>extxi[l+1]) break;
    }
  }
  for(int l=0;l<nbin-1;l++){
    if(yi[l]>0){
      xi[l]/=yi[l];
    }else{
      xi[l] = 0.5*(extxi[l] + extxi[l+1]);
    }
    yi[l] = yi[l]/dx/float(n);
  }
}
