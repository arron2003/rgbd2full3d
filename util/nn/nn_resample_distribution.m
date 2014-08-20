% Resamples the dataset (and labels) according to the given resampling
% strategy.
%
% Args:
%   nn
%   data - NxD
%   labels - Nx1
%
% Returns:
%   data - 
%   labels - 
function [data, labels, sampleWeights] = nn_resample_distribution(nn, data, labels, sampleWeights, minInsts)
  if nargin < 5
    minInsts = 10;
  end

  nn_consts;
  
  switch nn.resampleStrategy
    case RESAMPLE_RANDOM
      
      uniqueLabels = unique(labels);
      
      C = numel(uniqueLabels);
      minC = min(uniqueLabels);
      maxC = max(uniqueLabels);

      numInstances = histc(labels, [minC:maxC inf]);
      assert(numInstances(end) == 0);
      numInstances = numInstances(1:end-1);
      
      minInstances = min(numInstances(numInstances > minInsts));

      ndxs = [];

      for cc = 1 : C
        ndxsClass = find(labels == uniqueLabels(cc));

        % Take a random selection of these ndxs.
        seq = randperm(numel(ndxsClass));
        numSelect = floor(min([nn.resampleRatio * minInstances, numel(ndxsClass)]));

        ndxs = [ndxs; ndxsClass(seq(1:numSelect))];
      end
      
      data = data(ndxs, :);
      labels = labels(ndxs);
      sampleWeights = sampleWeights(ndxs);
      
    case RESAMPLE_BOOSTING
      
      % Select the 'most wrong' results.
      uniqueLabels = unique(labels);
      
      C = numel(uniqueLabels);
      minC = min(uniqueLabels);
      maxC = max(uniqueLabels);

      numInstances = histc(labels, [minC maxC inf]);
      assert(numInstances(end) == 0);
      numInstances = numInstances(1:end-1);
      
      minInstances = min(numInstances);

      % Now measure the error on each data point.
      outputs = nn_feed_forward(data, nn.weights, nn.biases, nn.outLayerType);
      errors = nn_calc_errors(outputs, labels, nn.outLayerType);
      
      ndxs = [];

      for cc = 1 : C
        ndxsClass = find(labels == uniqueLabels(cc));

        % Take the most mistaken instances of this label type.
        numToSelect = floor(min([nn.resampleRatio * minInstances, numel(ndxsClass)]));
        
        errorsForClass = errors(labels == uniqueLabels(cc));
        [~, seq] = sort(errorsForClass, 'descend');

        ndxs = [ndxs; ndxsClass(seq(1:numToSelect))];
      end
      
      data = data(ndxs, :);
      labels = labels(ndxs);
      sampleWeights = sampleWeights(ndxs);
  end
end
