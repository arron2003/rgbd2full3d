#include "mex.h"

#define ARG_REGIONS_TRUE 0
#define ARG_REGIONS_PRED 1
#define ARG_NUM_REGIONS_TRUE 2
#define ARG_NUM_REGIONS_PRED 3

void mexFunction(int nlhs, mxArray* plhs[],
                 const int nrhs, const mxArray* prhs[]) {
  
  if (nrhs != 4) {
    mexErrMsgTxt("4 arguments required.");
  }
  
  const int H = mxGetM(prhs[ARG_REGIONS_TRUE]);
  const int W = mxGetN(prhs[ARG_REGIONS_TRUE]);
  const int N = H * W;
  
  // First, double check the size of the images are the same.
  if (mxGetM(prhs[ARG_REGIONS_TRUE]) != mxGetM(prhs[ARG_REGIONS_TRUE]) ||
      mxGetN(prhs[ARG_REGIONS_TRUE]) != mxGetN(prhs[ARG_REGIONS_PRED])) {
    mexErrMsgTxt("imgRegionsTrue and imgRegionsPred must be the same size");
  }
  
  // Next, check that the data types are correct.
  if (mxGetClassID(prhs[ARG_REGIONS_TRUE]) != mxINT32_CLASS) {
    mexErrMsgTxt("imgRegionsTrue must be a 'int32'");
  } else if (mxGetClassID(prhs[ARG_REGIONS_PRED]) != mxINT32_CLASS) {
    mexErrMsgTxt("imgRegionsPRED must be a 'int32'");
  }

  int* regions_true = (int*) mxGetData(prhs[ARG_REGIONS_TRUE]);
  int* regions_pred = (int*) mxGetData(prhs[ARG_REGIONS_PRED]);
  
  int num_true_regions = mxGetScalar(prhs[ARG_NUM_REGIONS_TRUE]);
  int num_pred_regions = mxGetScalar(prhs[ARG_NUM_REGIONS_PRED]);

  plhs[0] = mxCreateDoubleMatrix(num_true_regions, num_pred_regions, mxREAL);
  double* overlap = (double*) mxGetData(plhs[0]);
  
  for (int ii = 0; ii < N; ++ii, ++regions_true, ++regions_pred) {
    int true_region_id = *regions_true;
    int pred_region_id = *regions_pred;
    
    if (true_region_id == 0 || pred_region_id == 0) {
      continue;
    }
    
    if (true_region_id > num_true_regions) {
      mexErrMsgTxt("True Region ID exceeds number of total regions.");
    } else if (pred_region_id > num_pred_regions) {
      mexErrMsgTxt("Pred Region ID exceeds number of total regions.");
    }    
    --true_region_id;
    --pred_region_id;
    
    ++overlap[true_region_id + pred_region_id * num_true_regions];
  }
}