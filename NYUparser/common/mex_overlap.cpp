#include "mex.h"

#define ARG_REGIONS_TRUE 0
#define ARG_REGIONS_PRED 1
#define ARG_REGION_IDS_TRUE 2
#define ARG_REGION_IDS_PRED 3

void mexFunction(int nlhs, mxArray* plhs[],
                 const int nrhs, const mxArray* prhs[]) {
  
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
  } else if (mxGetClassID(prhs[ARG_REGION_IDS_TRUE]) != mxINT32_CLASS) {
    mexErrMsgTxt("regionIdsTrue must be a 'int32'");
  } else if (mxGetClassID(prhs[ARG_REGION_IDS_PRED]) != mxINT32_CLASS) {
    mexErrMsgTxt("regionIdsPred must be a 'int32'");
  }

  int* regions_true = (int*) mxGetData(prhs[ARG_REGIONS_TRUE]);
  int* regions_pred = (int*) mxGetData(prhs[ARG_REGIONS_PRED]);
  
  int* region_ids_true = (int*) mxGetData(prhs[ARG_REGION_IDS_TRUE]);
  int* region_ids_pred = (int*) mxGetData(prhs[ARG_REGION_IDS_PRED]);
  
  int num_true_regions = static_cast<int>(mxGetScalar(prhs[ARG_REGION_IDS_TRUE]));
  int num_pred_regions = static_cast<int>(mxGetScalar(prhs[ARG_REGION_IDS_PRED]));

  const mwSize ndim = 2;
  mwSize dims[] = {num_true_regions, num_pred_regions};
  plhs[0] = mxCreateNumericArray(ndim, &dims[0], mxDOUBLE_CLASS, mxREAL);
  double* overlap = (double*) mxGetData(plhs[0]);
  
  plhs[1] = mxCreateNumericArray(ndim, &dims[0], mxDOUBLE_CLASS, mxREAL);
  double* cnt_isect = (double*) mxGetData(plhs[1]);
  
  plhs[2] = mxCreateNumericArray(ndim, &dims[0], mxDOUBLE_CLASS, mxREAL);
  double* cnt_union = (double*) mxGetData(plhs[2]);

  // Initialize both to 0.
  for (int ii = 0; ii < num_true_regions * num_pred_regions; ++ii) {
    cnt_union[ii] = 0;
    cnt_isect[ii] = 0;
  }
   
  // For each pixel.
  for (int nn = 0; nn < N; ++nn, ++regions_true, ++regions_pred) {
    int region_id_true = *regions_true;
    int region_id_pred = *regions_pred;

    if (region_id_true > 0 && region_id_pred > 0) {
      ++cnt_isect[ (region_id_pred-1) * num_true_regions + (region_id_true-1)];
    }

    if (region_id_pred > 0) {
      for (int ii = 0; ii < num_true_regions; ++ii) {
        ++cnt_union[(region_id_pred-1) * num_true_regions + ii];
      }
    }
    
    if (region_id_true > 0) {
      for (int jj = 0; jj < num_pred_regions; ++jj) {
        ++cnt_union[jj * num_true_regions + (region_id_true-1)];
      }
    }
  }
  
  for (int ii = 0; ii < num_true_regions * num_pred_regions; ++ii) {
    cnt_union[ii] -= cnt_isect[ii];
    
    if (cnt_union[ii] == 0) {
      overlap[ii] = 0;
    } else {
      overlap[ii] = static_cast<double>(cnt_isect[ii]) / cnt_union[ii];
    }
  }
}
