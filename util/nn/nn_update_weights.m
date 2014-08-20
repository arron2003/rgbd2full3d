function nn = nn_update_weights(nn, weightDerivs, biasDerivs)
  nn_consts;

  for ii = 1 : numel(nn.weights) - 1
    nn.weights{ii} = nn.weights{ii} - nn.eta * weightDerivs{ii} - nn.lambda * nn.weights{ii};
    nn.biases{ii} = nn.biases{ii} - nn.eta * biasDerivs{ii} - nn.lambda * nn.biases{ii};
  end
  
  % Now, for the final (top layer):
  if isfield(nn, 'classLambdas')
    assert(nn.outLayerType == SOFTMAX, ...
        'Class-based regularization is only valid with a SOFTMAX output');
    
    C = size(nn.weights{end}, 1);
      
    nn.weights{end} = nn.weights{end} - nn.eta * weightDerivs{end} - repmat(nn.classLambdas', [C 1]) .* nn.weights{end};
    nn.biases{end} = nn.biases{end} - nn.eta * biasDerivs{end} - nn.classLambdas' .* nn.biases{end};
  else
    nn.weights{end} = nn.weights{end} - nn.eta * weightDerivs{end} - nn.lambda * nn.weights{end};
    nn.biases{end} = nn.biases{end} - nn.eta * biasDerivs{end} - nn.lambda * nn.biases{end};
  end
end
