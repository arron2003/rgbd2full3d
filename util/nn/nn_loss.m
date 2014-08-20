% Calculates the loss for the given label/output pair and a specified loss
% function. 
%
% Args:
%   labels - NxC or Nx1 (depending on the lossFunc) matrix of labels.
%   outputs - NxC matrix
%   sampleWeights - Nx1 vector of weights for the samples.
%   lossFunc - the loss function ID. See nn_consts for the values.
%
% Returns:
%   totalLoss - the sum of the individual sample losses
%   sampleLoss - the loss (or cost) for each sample.
function [totalLoss, sampleLoss] = nn_loss(labels, outputs, sampleWeights, lossFunc)
  nn_consts;
  
  assert(~any(isnan(outputs(:))));

  if lossFunc == LOG_LOSS
    assert(ndims(labels) == 2);
    assert(size(labels, 2) == 1);
    
    sampleLoss = sampleWeights .* sum(-labels .* log(outputs+eps) - (1-labels) .* log(1-outputs+eps), 2);
  elseif lossFunc == MSE_LOSS
    assert(ndims(labels) == 2);
    
    sampleLoss = 0.5 * sampleWeights .* sum((labels - outputs) .^ 2, 2);
  elseif lossFunc == CROSS_ENTROPY_LOSS || lossFunc == CROSS_ENTROPY_RANK_LOSS
    % 2 classes
    %loss = -sum(sum(labels .* log(outputs) + (1-labels) .* log(1-outputs)));

    % multiclass
    sampleLoss = -sampleWeights .* sum(labels .* log(eps + outputs ./ (labels+eps)), 2);
    
    if isa(sampleLoss, 'parallel.gpu.GPUArray')
      assert(gather(sum(isnan(sampleLoss(:)))) == 0);
    else
      assert(sum(isnan(sampleLoss(:))) == 0);
    end
  end
  
  totalLoss = sum(sampleLoss);
end
