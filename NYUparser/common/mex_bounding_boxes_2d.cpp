#include "mex.h"

#define ARG_REGIONS 0
#define ARG_REGION_IDS 1

void mexFunction(int nlhs, mxArray* plhs[],
                 const int nrhs, const mxArray* prhs[]) {
  
  if (nrhs != 2) {
    mexErrMsgTxt("Exactly 2 arguments required.");
  }
  
  const int H = mxGetM(prhs[ARG_REGIONS]);
  const int W = mxGetN(prhs[ARG_REGIONS]);
  const int N = H * W;
  
  // Next, check that the data types are correct.
  if (mxGetClassID(prhs[ARG_REGIONS]) != mxINT32_CLASS) {
    mexErrMsgTxt("imgRegion must be a 'int32'");
  } else if (mxGetClassID(prhs[ARG_REGION_IDS]) != mxINT32_CLASS) {
    mexErrMsgTxt("regionIds must be a 'int32'");
  }
  
  int* regions = (int*) mxGetData(prhs[ARG_REGIONS]);
  int R = mxGetNumberOfElements(prhs[ARG_REGION_IDS]);
  
  // X min, X max, Y min, Y max
  plhs[0] = mxCreateDoubleMatrix(4, R, mxREAL);
  double* bb = (double*) mxGetData(plhs[0]);
  
  for (int ii = 0; ii < R; ++ii) {
    bb[4 * ii + 0] = 10e10;
    bb[4 * ii + 1] = 0;
    bb[4 * ii + 2] = 10e10;
    bb[4 * ii + 3] = 0;
  }
  
  for (int ii = 0; ii < N; ++ii, ++regions) {
    int region_id = *regions;
    if (region_id == 0) {
      continue;
    }
    
    if (region_id > R) {
      mexErrMsgTxt("Region ID exceeds number of total regions.");
    } 
    --region_id;
    
    // Adding 1 for matlab 1-indexing.
    int Y = (ii % H) + 1;
    int X = (ii / H) + 1;
    
    if (bb[4 * region_id + 0] > X) {
      bb[4 * region_id + 0] = X;
    }
    if (bb[4 * region_id + 1] < X) {
      bb[4 * region_id + 1] = X;
    }
    
    if (bb[4 * region_id + 2] > Y) {
      bb[4 * region_id + 2] = Y;
    }
    if (bb[4 * region_id + 3] < Y) {
      bb[4 * region_id + 3] = Y;
    }
  }
}