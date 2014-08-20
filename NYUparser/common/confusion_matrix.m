% Calculates the confusion given a set of labels and predictions.
%
% Args:
%   labels - the set of labels which range from 1..C, an Nx1 vector.
%   predictions - the set of predictions which range from 1..C, an 
%                 Nx1 vector.
%   numClasses - the maximum number of classes.
function confMat = confusion_matrix(labels, predictions, numClasses)
  error(nargchk(3,3,nargin));
  labels = uint16(labels(:));
  predictions = uint16(predictions(:));
  numClasses = uint16(numClasses);
  confMat = mex_conf_mat(labels, predictions, numClasses);
end
