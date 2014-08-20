% Evaluates the accuracy of a given segmentation using pixel-level
% accuracy.
%
% Args:
%   predictions - a HxW matrix of labels, where H and W are the height and
%                 width of the image, respectively.
%   groundTruth - a HxW matrix of labels, where H and W are the height and
%                 width of the image, respectively.
%   numClasses - the number of total classes in the experiment.
function [accuracy, numCorrect, numLabeled, cm] = eval_seg(predictions, ...
    groundTruth, numClasses)
  
  % Predictions should only be 0 in the case where we are evaluating
  % Ground-Truth-region labels (where some of the pixels are not labeled in
  % the ground-truth case).
  labeledNdxs = groundTruth > 0 & predictions > 0;
  numLabeled = nnz(labeledNdxs);
  
  groundTruth = groundTruth(labeledNdxs);
  predictions = predictions(labeledNdxs);
  
  numCorrect = nnz(predictions == groundTruth);
  
  accuracy = numCorrect / numel(groundTruth);
  
  cm = confusion_matrix(groundTruth, predictions, numClasses);
end
