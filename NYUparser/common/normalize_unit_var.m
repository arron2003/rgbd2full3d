% Normalizes the data so that each dimension has unit variance.
%
% Args:
%   data - NxD matrix
% 
% Returns:
%   data - NxD matrix
%   stds - 1xD matrix of standard deviations prior to normalization.
function [data, stds] = normalize_unit_var(data, stds)
  N = size(data, 1);
  if nargin == 1
    stds = std(data, [], 1);
  end
  data = data ./ (eps + ones(N, 1) * stds);
end