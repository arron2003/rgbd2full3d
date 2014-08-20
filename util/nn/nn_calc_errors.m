% Calculates the errors to propagate back from the top layer, e.g:
%  dE/dy where E is the loss function and y represents the output layer.
% 
% Args:
%   outputs - NxD matrix where N is the number of samples and D is the
%             dimensionaliaty of the output.
%   labels - NxD matrix where N is the number of samples and D is the
%            dimensionality of the output.
%   outLayerType - the type of output layer being used. See rbm_consts.
%
% Returns:
%   errors - NxD matrix where N is the number of samples and D is the
%            dimensionaliaty of the output.
function errors = nn_calc_errors(outputs, labels, outLayerType)
  nn_consts;
  switch outLayerType
    case BINARY
      errors = -labels ./ (outputs+eps) + (1-labels) ./ (1-outputs + eps);
    case GAUSSIAN
      errors = outputs - labels;
    case SOFTMAX
      errors = outputs - labels;
    otherwise
      error('OutLayerType %d not recognized.', outLayerType);
  end
end