% Evaluates the prediction of a given neural network.
%
% Args:
%   data - NxD matrix where N is the number of data points and D is the
%          dimensionality of the data.
%   labels - NxC matrix of 'verbose' labels.
%
% Returns:
%   accuracy - the number of correct matches divided by total samples.
%   confMat - confusion matrix.
%   ranks - the ranks of the correct labels.
%   outputs - the outputs from feeding forward
function [accuracy, confMat, ranks, outputs] = nn_eval(nn, data, labels)
  nn_consts;
  
  ranks = [];
  
  C = numel(nn.biases{end});
  
  if nn.outLayerType == SOFTMAX
    targets = labels;
    labels = get_verbose_labels(labels, 1:C);
  end

  outputs = nn_feed_forward(nn, data);
  
  % Reshape the outputs and labels if they are batched.
  [N, D, B] = size(outputs);
  assert(size(labels, 1) == N);
  assert(size(labels, 3) == B);
  
  if size(labels, 2) > 1
    C = size(labels, 2);
    labels = reshape(permute(labels, [1 3 2]), [N*B C]);
  else
    C = numel(unique(labels));
  end
  
  outputs = reshape(permute(outputs, [1 3 2]), [N*B D]);
  
  if nn.outLayerType == SOFTMAX
    [~, inds] = sort(outputs, 2, 'descend');
    [~, predictions] = max(outputs, [], 2);
    
    ranks = zeros(N*B, 1);
    
    for ii = 1 : N*B
      ranks(ii) = find(inds(ii, :) == targets(ii));
    end
  elseif nn.outLayerType == BINARY || nn.outLayerType == GAUSSIAN
    predictions = (outputs >= 0.5);
  end

  if size(labels,2) > 1
    [~, labels] = max(labels, [], 2);
  end
  labels = labels(:);

  accuracy = sum(predictions == labels) / numel(labels);
  
  if nn.outLayerType == BINARY
    labels = labels + 1;
    predictions = predictions + 1;
  end
  
  confMat = confusion_matrix(labels, predictions, C);
end
