function [verboseLabels, targets] = get_verbose_labels(simpleLabels, alphabet)
% Generates a set of binary labels, each of which is an L-vector where
% L is the cardinality of the set of labels. The i-th entry of the L-vector
% representing the n-th simple label is equal to 1 if simpleLabels(n) == i
% and 0 otherwise.
%
% For example, given the simple labels [0 1 1 2 0], the label alphabet
% consists of [0 1 2], and the verbose labels will be:
%
% [1 0 0;
%  0 1 0;
%  0 1 0;
%  0 0 1;
%  1 0 0]
%
% Inputs:
%   simpleLabels - a matrix of discrete labels in any order.
%   alphabet - an optional argument. If the alphabet is not provided, one
%              is inferred based on the unique values contained within
%              simpleLabels.
%
% Outputs:
%   verboseLabels - a set of verboseLabels representing the simpleLabels
%                   given an alphabet.
%   targets - a new set of simple labels. This will be different from the
%             simpleLabels provided if the simpleLabels do not contain
%             integer labels 1..N. For example, if the simpleLabels contain
%             [0 0.5 0.75], the targets will be [1 2 3].

error(nargchk(1, 2, nargin));

if nargin == 1
  alphabet = unique(simpleLabels(:));
end

originalSize = size(simpleLabels);
alphabet = alphabet(:);
simpleLabels = simpleLabels(:);

verboseLabels = zeros(length(alphabet), length(simpleLabels));

targets = zeros(1, numel(simpleLabels));

for i = 1 : length(alphabet)
  targets(simpleLabels == alphabet(i)) = i;
  verboseLabels(i, targets == i) = 1;
end

targets = reshape(targets, originalSize);
verboseLabels = squeeze(reshape(verboseLabels, [length(alphabet) originalSize]));
verboseLabels = verboseLabels';
