// Returns masks for each regions' border (edge) and interior components.
#include "mex.h"

#define ARG_IMG_REGIONS 0
#define ARG_NUM_REGIONS 1
#define ARG_IMG_ERODED 2

#define RET_INTERIORS 0
#define RET_EDGES 1

void mexFunction(int nlhs, mxArray* plhs[],
                 const int nrhs, const mxArray* prhs[]) {
  
  if (nrhs != 3) {
    mexPrintf("Exactly 3 arguments required.");
  }
  
  // Grab the region mask (HxW).
  int H = mxGetM(prhs[ARG_IMG_REGIONS]);
  int W = mxGetN(prhs[ARG_IMG_REGIONS]);
  int N = H * W;
  
  if (mxGetClassID(prhs[ARG_IMG_REGIONS]) != mxINT32_CLASS) {
    mexErrMsgTxt("imgRegions must be of class 'int32'");
  }
  int* img_regions = (int*) mxGetData(prhs[ARG_IMG_REGIONS]);
  
  // Get the number of regions.
  int R = (int) mxGetScalar(prhs[ARG_NUM_REGIONS]);
  
  // Grab the eroded region mask.
  if (mxGetNumberOfDimensions(prhs[ARG_IMG_ERODED]) != 2
      || mxGetNumberOfElements(prhs[ARG_IMG_ERODED]) != N) {
    mexErrMsgTxt("imgEroded must be of size HxW");
  } else if (mxGetClassID(prhs[ARG_IMG_ERODED]) != mxLOGICAL_CLASS) {
    mexErrMsgTxt("imgEroded must be of class 'logical'");
  }
  bool* img_eroded = (bool*) mxGetData(prhs[ARG_IMG_ERODED]);
  
  // Create the outputs.
  mwSize ndims_out = 3;
  mwSize dims_out[] = {H, W, R};
  plhs[RET_INTERIORS] = mxCreateLogicalArray(ndims_out, &dims_out[0]);
  plhs[RET_EDGES] = mxCreateLogicalArray(ndims_out, &dims_out[0]);
  
  mwSize ndims_out_size = 2;
  mwSize dims_out_size[] = {R, 1};

  bool* img_interiors = (bool*) mxGetData(plhs[RET_INTERIORS]);
  bool* img_edges = (bool*) mxGetData(plhs[RET_EDGES]);
  
  unsigned int num_pix_interiors[R];
  unsigned int num_pix_edges[R];
  
  for (int rr = 0; rr < R; ++rr) {
    num_pix_interiors[rr] = 0;
    num_pix_edges[rr] = 0;
  }

  for (int ii = 0; ii < N; ++ii, ++img_eroded) {
    // Grab the Region ID.
    int region_id = img_regions[ii] - 1;
    
    // Is it in the interior?
    if (*img_eroded) {
      img_interiors[region_id * N + ii] = true;
      ++num_pix_interiors[region_id];
    } else {
      img_edges[region_id * N + ii] = true;
      ++num_pix_edges[region_id];
    }
  }
  
  // Now that we know how many pixels make up the interiors and borders,
  // if we have any interiors or border regions that are empty, just use
  // the entire mask instead.
  for (int ii = 0; ii < N; ++ii) {
    // Grab the Region ID.
    int region_id = img_regions[ii] - 1;
    
    if (num_pix_interiors[region_id] == 0) {
      img_interiors[region_id * N + ii] = true;
    }
    
    if (num_pix_edges[region_id] == 0) {
      img_edges[region_id * N + ii] = true;
    }
  }
}
