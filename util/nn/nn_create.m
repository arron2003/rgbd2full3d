% Creates the neural network.
%
% Args:
%   nn - the neural network struct.
%   tag - (optional) a name for the neural network.
%
% Returns:
%   nn - the neural network.
function nn = nn_create(dims, tag)
  nn_consts;

  nn = struct();
  nn.debug = 0;
  nn.dims = dims;
  nn.lambda = 0;
  nn.batchSize = 100;
  
  nn.actFunc = LOGISTIC;
  nn.resampleStrategy = RESAMPLE_NONE;
  nn.resampleRatio = 0;
  
  nn.rankSig = 10;
  
  if nargin < 2
    nn.tag = '';
  else
    nn.tag = tag;
  end
  
  nn.weights = cell(0);
  nn.biases = cell(0);
end
