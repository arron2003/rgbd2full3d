% Normalizes the data so that each dimension has zero mean.
%
% Args:
%   data - NxD matrix
% 
% Returns:
%   data - NxD matrix
%   means - 1xD matrix of means prior to normalization.
function [data, means] = normalize_zero_mean(data, means)
  N = size(data, 1);
  if nargin == 1
    means = mean(data, 1);
  end
  data = data - ones(N, 1) * means;
end