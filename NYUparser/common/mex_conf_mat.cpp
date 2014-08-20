// Creates a confusion matrix given a series of labels and predictions.
// Note that the code assumes that the labels range from 1..N (1-indexed, 
// and not 0-indexed).

#include "mex.h"

#define ARG_LABELS 0
#define ARG_PREDS 1
#define ARG_CLASSES 2

// Throws an exception if the given vector is invalid.
void validate_vec(const mxArray* vec) {
  if (mxGetNumberOfDimensions(vec) != 2 || mxGetN(vec) > 1) {
    mexErrMsgTxt("labels and predictions must be a column vector");
  }
  
  if (mxGetClassID(vec) != mxUINT16_CLASS) {
    mexErrMsgTxt("labels and predictions must be uint16");
  }
}

// Assumes that the inputs range from 1 to N.
void mexFunction(int nlhs, mxArray* plhs[],
                 const int nrhs, const mxArray* prhs[]) {  
  if (nrhs != 3) {
    mexErrMsgTxt("Usage: conf_mat(labels, predictions, num_classes)");
  }
  
  validate_vec(prhs[ARG_LABELS]);
  validate_vec(prhs[ARG_PREDS]);
  
  // Make sure they have the same length.
  if (mxGetNumberOfElements(prhs[ARG_LABELS]) != 
      mxGetNumberOfElements(prhs[ARG_PREDS])) {
    mexErrMsgTxt("Number of labels and predictions MUST match");
  }
  
  const int N = mxGetNumberOfElements(prhs[ARG_LABELS]);
  const int num_classes = (int) mxGetScalar(prhs[ARG_CLASSES]);
  if (num_classes <= 1) {
    mexErrMsgTxt("Number of classes must exceed 1");
  }
  
  // Allocate the output.
  plhs[0] = mxCreateDoubleMatrix(num_classes, num_classes, mxREAL);
  double* conf_mat = (double*) mxGetData(plhs[0]);
  
  unsigned short int* labels = (unsigned short int*) mxGetData(prhs[ARG_LABELS]);
  unsigned short int* preds = (unsigned short int*) mxGetData(prhs[ARG_PREDS]);
  
  for (int ii = 0; ii < N; ++ii, ++labels, ++preds) {
    unsigned short int lbl = *labels;
    unsigned short int pred = *preds;
    
    if (lbl == 0 || lbl > num_classes) {
      mexErrMsgIdAndTxt("params:labels", "Label(%d) is out of range", ii + 1);
    } else if (pred == 0 || pred > num_classes) {
      mexErrMsgIdAndTxt("params:predictions", "Prediction(%d) is out of range", ii + 1);
    }
    
    --lbl;
    --pred;
    
    int abs_ndx = pred * num_classes + lbl;
    ++conf_mat[abs_ndx];
  }
}
