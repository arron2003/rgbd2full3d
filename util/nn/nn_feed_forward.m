% Evaluates the neural network by feeding forward the given data.
%
% Args:
%   nn - the neural network.
%   data - a NxD matrix where N is the number of data points and D is the
%          dimensionality of the data
%
% Returns:
%   output - the outer activations of the network
%   units - the activations of every layer of the network.
function [outputs, acts] = nn_feed_forward(nn, data)
  nn_consts;
  N = size(data, 1);
  
  assert(numel(nn.weights) == numel(nn.biases));
  
  numLayers = numel(nn.weights);
  
  acts = cell(numLayers, 1);
  
  for ii = 1 : numLayers
    if ii == 1
      prevActs = data;
    else
      prevActs = acts{ii-1};
    end

    acts{ii} = prevActs * nn.weights{ii} + ones(N, 1) * nn.biases{ii};

    if ii == numLayers && nn.outLayerType == SOFTMAX
      % Ugly hack to fix inf issues. 
      % TODO: figure out why this is happening.
      acts{ii}(acts{ii} > 80) = 80;
      
      acts{ii} = exp(acts{ii});
      denom = repmat(sum(acts{ii}, 2), 1, size(acts{ii}, 2));
      acts{ii} = acts{ii} ./ (denom+eps);
    elseif ii == numLayers && nn.outLayerType == GAUSSIAN
      % do nothing in this case.
    else
      acts{ii} = 1 ./ (1 + exp(-acts{ii}));
    end
  end

  outputs = acts{end};
end
