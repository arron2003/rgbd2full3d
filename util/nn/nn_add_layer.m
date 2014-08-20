% Adds a layer to the neural network.
%
% Args:
%   nn - the neural network struct.
%   numHid - the number of hidden units in the layer being added.
%
% Returns: 
%   nn - the neural network, with the added layer.
function nn = nn_add_layer(nn, numHid, outLayerType, lossFunc, weightCoeff)
  nn_consts;
  assert(isstruct(nn));
  assert(isfield(nn, 'dims'), 'Whoops, nn.dims is missing');
  assert(isscalar(numHid));
  
  if nargin < 3
    nn.outLayerType = BINARY;
  else
    nn.outLayerType = outLayerType;
  end
  
  if nargin < 4
    switch nn.outLayerType
      case BINARY
        nn.lossFunc = LOG_LOSS;
      case GAUSSIAN
        nn.lossFunc = MSE_LOSS;
      case SOFTMAX
        nn.lossFunc = CROSS_ENTROPY_LOSS;
    end
  else
    nn.lossFunc = lossFunc;
  end
  
  if nargin < 5
    weightCoeff = .01;
  end
  
  if isfield(nn, 'biases') && ~isempty(nn.biases);
    L = numel(nn.weights);
    dims = numel(nn.biases{end});
  else
    L = 0;
    dims = nn.dims;
  end
  
  nn.weights{L+1} = weightCoeff * randn(dims, numHid);
  nn.biases{L+1} = weightCoeff * randn(1, numHid);
end
