% Performs backpropagating over a given network.
%
% Args:
%   errors - the error term at the top of the network.
%   units - an array of cells, where each cell contains the activations of a
%           given layer. Each cell should contain a matrix of size NxD
%           where N is then number of samples (think number of inputs to
%           the network) and D is the size of the layer, or the input size
%           for units{1}.
%   weights - an array of weights
%   outputLayerType - the type of output layer (see nn_consts.m);
function [df dw db] = ...
    nn_backpropagate(errors, inputData, units, weights, outputLayerType, sampleWeights)
 
  nn_consts;
  if nargin < 6
    sampleWeights = ones(size(errors, 1), 1);
  end
  
  % The errors that are backpropagated back. More formally, this is the
  % dE/Da_j term where E is the error function and a_j is the weighted
  % input to unit a_j.
  d = {};
  
  % The final partial derivative of the weights. More formally, this is the
  % dE/Dw_ij term where E is the error function and w_ij is the weight
  % connecting units i and j.
  dw = {};
  
  % The final partial derivative of the biases. More formally, this is the
  % dE/Dw_j term where E is the error function and w_j is the bias on unit
  % j.
  db = {};
  
  for ii = length(units) : -1 : 1
    if ii == 1
      prevActivations = inputData;
    else
      prevActivations = units{ii-1};
    end

    if ii == length(units) && outputLayerType == BINARY
      d{ii} = units{ii} .* (1-units{ii}) .* errors;
      dw{ii} = prevActivations' * bsxfun(@times, d{ii}, sampleWeights);
      db{ii} = sum(bsxfun(@times, d{ii}, sampleWeights), 1);
    elseif ii == length(units) && outputLayerType == SOFTMAX
      d{ii} = errors;
      dw{ii} = prevActivations' * bsxfun(@times, d{ii}, sampleWeights);
      db{ii} = sum(bsxfun(@times, d{ii}, sampleWeights), 1);

      % Use this code if you use the euclidean loss function instead of
      % cross entropy.
%       numhid = size(errors,2);
%       d{ii} = (errors - repmat(diag(errors*units{ii}'),1,numhid)).*units{ii};
    elseif ii == length(units) && outputLayerType == GAUSSIAN
      d{ii} = errors;
      dw{ii} = prevActivations' * bsxfun(@times, d{ii}, sampleWeights);
      db{ii} = sum(bsxfun(@times, d{ii}, sampleWeights), 1);
    else
      d{ii} = (units{ii} .* (1-units{ii})) .* (weights{ii+1} * d{ii+1}')';
      dw{ii} = prevActivations' * bsxfun(@times, d{ii}, sampleWeights);
      db{ii} = sum(bsxfun(@times, d{ii}, sampleWeights), 1);
    end
  end

  df = nn_vectorize_weights(dw, db);
end