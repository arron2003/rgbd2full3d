function [f, df, outputs, errors, sampleLoss] = nn_grad_descent(weightVector, nn, data, labels, sampleWeights)
  nn_consts;

  if nargin < 5
    sampleWeights = ones(size(data,1),1);
  end
  
  
  [nn.weights, nn.biases] = nn_devectorize_weights(nn, weightVector);
  
  [outputs, units] = nn_feed_forward(nn, data);
  errors = nn_calc_errors(outputs, labels, nn.outLayerType);
  
  
  % Update the sample weights
  switch nn.lossFunc
    case CROSS_ENTROPY_RANK_LOSS
      rankLossWeight = get_rank_aware_weights(labels, outputs, nn.maxAllowableRank, nn.rankLossWeightType, nn.rankSig);
      sampleWeights = sampleWeights .* rankLossWeight;
  end
  
  df = nn_backpropagate(errors, data, units, nn.weights, nn.outLayerType, sampleWeights);
  
  [f, sampleLoss] = nn_loss(labels, outputs, sampleWeights, nn.lossFunc);
end

function rankLossWeight = get_rank_aware_weights(labels, outputs, allowableRank, rankLossWeightType, rankSig)
  nn_consts; 

  % Figure out in which samples the rank of the correct answer is "good
  % enough".

  % First, figure out the expected rank.
  [~, targets] = max(labels, [], 2);

  % Next, measure the observed rank.
  [~, inds] = sort(outputs, 2, 'descend');

  % Calculate parts of the loss function.
  switch rankLossWeightType
    case RANK_LOSS_LINEAR
     m = 1 / allowableRank;
     b = -m;

     weights = m * (1:allowableRank) + b;
     assert(all(weights <= 1));
    case RANK_LOSS_EXP
      weights = exp(-(5:-1:1)./rankSig);
      weights = weights ./ sum(weights);
    case RANK_LOSS_TRUNCATED
      weights = zeros(allowableRank, 1);
  end

  N = size(outputs, 1);
  rankLossWeight = ones(N, 1);
  for ii = 1 : N
    rank = find(targets(ii) == inds(ii, :));
    if rank > allowableRank
      continue;
    end

    switch rankLossWeightType
      case RANK_LOSS_TRUNCATED
        rankLossWeight(ii) = 0;
      case RANK_LOSS_LINEAR
        rankLossWeight(ii) = weights(rank);
      case RANK_LOSS_EXP
        rankLossWeight(ii) = weights(rank);
      otherwise
        error('not supported');
    end
  end
end
