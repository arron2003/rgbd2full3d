% Calculates and plots the ROC curve for the given predictions and labels.
%
% Args:
%   scores - Nx1 vector of confidences or scores.
%   labels - Nx1 vector of binary labels.
%
% Returns:
%   auc - the auc under the ROC curve
function auc = roc(scores, labels, color)
  assert(ndims(scores) == ndims(labels));
  assert(all(size(scores) == size(labels)));
  assert(islogical(labels), 'Labels must be logical');

  [~, inds] = sort(scores, 'descend');
  truePositives = cumsum(labels(inds));
  falsePositives = cumsum(~labels(inds));

  if nargin < 3
    color = '--b';
  else
    assert(numel(color) == 1);
    color = sprintf('--%s', color);
  end
  
  X = falsePositives / falsePositives(end);
  Y = truePositives / truePositives(end);
  plot(X, Y, color);
  
  axis([0 1 0 1]);
  grid on;
  
  % Calculate the auc under the curve.
  auc = 0;
  step = 100;
  for ii = 1 + step : step : numel(scores)
    auc = auc + (X(ii)-X(ii-step)) * (Y(ii)+Y(ii-step))/2;
  end
end
