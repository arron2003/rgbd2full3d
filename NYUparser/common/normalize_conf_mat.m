% Normalizes the given confusion matrix.
%
% Args:
%   cm - CxC confusion matrix.
% 
% Returns:
%   cm - the normalized confusion matrix such that each row sums to 1.
function cm = normalize_conf_mat(cm)
  assert(ndims(cm) == 2);
  assert(size(cm, 1) == size(cm, 2));
  assert(all(all(cm >= 0)));

  C = size(cm, 1);
  
  cm = cm ./ repmat(sum(cm, 2), [1 C]);
  cm(isnan(cm)) = 0;
end