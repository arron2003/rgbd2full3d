% Shuffles the given data set.
%
% Args:
%   data - NxD matrix where N is the number of samples and D is the
%          dimensionality of the data.
%
% Returns:
%   data - NxD matrix where the samples have been shuffled.
%   seq - the sequence used to randomize the data.
function [data, seq] = shuffle_data(data)
  N = size(data, 1);
  seq = randperm(N);
  data = data(seq, :);
end