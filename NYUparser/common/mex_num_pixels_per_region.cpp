#include "mex.h"

#define ARG_REGIONS 0
#define ARG_NUM_REGIONS 1

void mexFunction(int nlhs, mxArray* plhs[],
                 const int nrhs, const mxArray* prhs[]) {
  
  if (nrhs != 2) {
    mexErrMsgTxt("2 arguments required.");
  }
  
  const int H = mxGetM(prhs[ARG_REGIONS]);
  const int W = mxGetN(prhs[ARG_REGIONS]);
  const int N = H * W;
  
  // Next, check that the data types are correct.
  if (mxGetClassID(prhs[ARG_REGIONS]) != mxINT32_CLASS) {
    mexErrMsgTxt("imgRegionsTrue must be a 'int32'");
  }

  int* regions = (int*) mxGetData(prhs[ARG_REGIONS]);
  int R = (int) mxGetScalar(prhs[ARG_NUM_REGIONS]);

  plhs[0] = mxCreateDoubleMatrix(R, 1, mxREAL);
  double* pixel_counts = (double*) mxGetData(plhs[0]);
  
  for (int ii = 0; ii < N; ++ii, ++regions) {
    int region_id = *regions;
    
    if (region_id == 0) {
      continue;
    }
    
    if (region_id > R) {
      mexErrMsgTxt("Region ID exceeds number of total regions.");
    }
    
    --region_id;
    
    ++pixel_counts[region_id];
  }
}