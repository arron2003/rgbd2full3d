% The orderering is from bottom layer to top layer and all weights come
% before all biases.
% 
% Args:
%   weightsAndBiases - vector of weights and biases.
function [weights biases] = nn_devectorize_weights(nn, weightsAndBiases)
  assert(numel(weightsAndBiases) == length(weightsAndBiases));
  pointer = 1; % Keeps track of where we are linearly in the weight vector.
  
  weights = cell(numel(nn.weights),1);
  biases = cell(numel(nn.biases),1);
  
  for ii = 1 : length(nn.weights)
    weights{ii} = weightsAndBiases(pointer : pointer + numel(nn.weights{ii}) - 1);
    weights{ii} = reshape(weights{ii}, size(nn.weights{ii}));
    pointer = pointer + numel(nn.weights{ii});
  end
  for ii = 1 : length(nn.biases)
    biases{ii} = weightsAndBiases(pointer : pointer + numel(nn.biases{ii}) - 1);
    biases{ii} = reshape(biases{ii}, size(nn.biases{ii}));
    pointer = pointer + numel(nn.biases{ii});
  end
end