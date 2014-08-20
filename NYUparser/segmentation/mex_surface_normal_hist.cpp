#include "mex.h"

#define ARG_IMG_REGIONS 0
#define ARG_NUM_REGIONS 1
#define ARG_INCLINATIONS 2
#define ARG_INC_BINS 3
#define ARG_AZIMUTHS 4
#define ARG_AZ_BINS 5

#define RET_INC_HISTS 0
#define RET_AZ_HISTS 1

void mexFunction(int nlhs, mxArray* plhs[],
                 const int nrhs, const mxArray* prhs[]) {
  
  if (nrhs != 6) {
    mexErrMsgTxt("6 args required.");
  }
  
  int* img_regions = (int*) mxGetData(prhs[ARG_IMG_REGIONS]);
  int N = (int) mxGetScalar(prhs[ARG_NUM_REGIONS]);
  double* inclinations = (double*) mxGetData(prhs[ARG_INCLINATIONS]);
  double* bins_inc = (double*) mxGetData(prhs[ARG_INC_BINS]);
  double* azimuths = (double*) mxGetData(prhs[ARG_AZIMUTHS]);
  double* bins_az = (double*) mxGetData(prhs[ARG_AZ_BINS]);
  
  int num_pix = mxGetNumberOfElements(prhs[ARG_IMG_REGIONS]);
  int num_bins_inc = mxGetNumberOfElements(prhs[ARG_INC_BINS]);
  int num_bins_az = mxGetNumberOfElements(prhs[ARG_AZ_BINS]);
  
  plhs[RET_INC_HISTS] = mxCreateDoubleMatrix(num_bins_inc, N, mxREAL);
  plhs[RET_AZ_HISTS] = mxCreateDoubleMatrix(num_bins_az, N, mxREAL);
  
  double* inc_hists = (double*) mxGetData(plhs[RET_INC_HISTS]);
  double* az_hists = (double*) mxGetData(plhs[RET_AZ_HISTS]);
  
  
  for (int ii = 0; ii < num_pix; ++ii, ++img_regions, ++inclinations, ++azimuths) {
    int region_id = *img_regions;
    if (region_id == 0) {
      continue;
    }
    
    --region_id;
    
    // Figure out which inclination bin to fall into.
    double inclination = *inclinations;
    int inc_bin = -1;
    for (int jj = 0; jj < num_bins_inc; ++jj) {
      if (inclination >= bins_inc[jj]) {
        inc_bin = jj;
      }
    }
    
    ++inc_hists[region_id * num_bins_inc + inc_bin];
    
    // Figure out which azimuth bin to fall into.
    double azimuth = *azimuths;
    int az_bin = -1;
    for (int jj = 0; jj < num_bins_az; ++jj) {
      if (azimuth >= bins_az[jj]) {
        az_bin = jj;
      }
    }
    ++az_hists[region_id * num_bins_az + az_bin];
  }  
}