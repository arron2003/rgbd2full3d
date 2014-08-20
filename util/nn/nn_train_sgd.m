% Trains the NeuralNetwork
%
% Inputs:
%   data - a NxD matrix where D is the dimensionality of the training data,
%          and N is the number of samples.
%   labels - a NxC matrix where N is the number of samples and C is the
%            number of classes.
%   sampleWeights - (optional) Nx1 matrix used to weight the backpropagated
%                   derivatives. 
function nn = nn_train_sgd(nn, data, labels, sampleWeights)
  nn_consts;
  assert(size(data, 1) == size(labels, 1));
  
  if nargin < 4
    sampleWeights = ones(size(data, 1), 1);
  end
  
  C = numel(nn.biases{end});
  
  if ~isfield(nn, 'updatesMade')
    nn.updatesMade = 0;
  end
  
  % Get the Mini Batch size.
  if ~isfield(nn, 'batchSize')
    fprintf('Mini-batch size not set. Using default of 100.');
    batchSize = 100;
  else
    batchSize = nn.batchSize;
  end
  
  updatesToMake = nn.numUpdates - nn.updatesMade;
  
  epoch = 0;
  sse = zeros(0, 1);
  epochTimes = zeros(0, 1);
  startTimeRun = clock;
  
  while updatesToMake > 0
    epoch = epoch + 1;
    sse(epoch) = 0;
    epochTimes(epoch) = 0;
    
    startTimeEpoch = clock; % Starts counting the time of the current epoch.
    
    [dataSample, labelsSample, weightsSample] = nn_resample_distribution(nn, data, labels, sampleWeights);

    % Shuffle the dataset.
    [dataSample, seq] = shuffle_data(dataSample);
    labelsSample = labelsSample(seq,:);
    weightsSample = weightsSample(seq);
    
    % Make the labels verbose
    if nn.outLayerType == SOFTMAX
      labelsSample = get_verbose_labels(labelsSample, 1:C);
    end

    N = size(dataSample, 1);

    numBatches = ceil(N / batchSize);

    for bb = 1 : numBatches
      startNdx = (bb-1) * batchSize + 1;
      endNdx = min([bb * batchSize, N]);

      dataBatch = dataSample(startNdx:endNdx, :);
      labelBatch = labelsSample(startNdx:endNdx, :);
      weightsSampleBatch = weightsSample(startNdx:endNdx, :);
      
      numelBatch = endNdx - startNdx + 1;

      W = nn_vectorize_weights(nn.weights, nn.biases);
      [loss, df, outputs] = nn_grad_descent(W, nn, ...
          dataBatch, labelBatch, weightsSampleBatch);

      assert(~isnan(loss));
        
      avgLoss = loss / numelBatch;
      sse(epoch) = sse(epoch) + loss;

      % fprintf('(%s) Batch %d/%d, avgLoss=%f\r', nn.tag, bb, numBatches, avgLoss); 
      [weightDerivs biasDerivs] = nn_devectorize_weights(nn, df);

      nn = nn_update_weights(nn, weightDerivs, biasDerivs);

      nn.updatesMade = nn.updatesMade + 1;
      if nn.updatesMade >= nn.numUpdates
        break;
      end
    end

    mseEpoch = sse(epoch) / endNdx;
    % fprintf('(%s) Updates %d/%d, SSE=%f, MSE: %f\n', nn.tag, nn.updatesMade, nn.numUpdates, sse(epoch), mseEpoch);

    if nn.debug == 1
      f = sfigure(1);
      subplot(3, 2, 1); imagesc(single(labelBatch)); colormap('gray'); title('labels');
      subplot(3, 2, 2); imagesc(single(outputs)); colormap('gray'); title('predictions'); caxis([0 1]);
      subplot(3, 2, 3); plot(max(1, epoch-500):epoch, sse(max(1, epoch-500):epoch)); title('Errors');
      subplot(3, 2, 4); imagesc(single(nn.weights{1})); title('1st layer weights');
      subplot(3, 2, 5); hist(nn.weights{1}(:), 20); title('Hist of Weights');
      subplot(3, 2, 6); hist(nn.biases{1}(:), 20); title('Hist of Biases');
      set(f, 'Name', sprintf('NN - %s', nn.tag));
      drawnow;
    end
    
    epochTimes(epoch) = etime(clock, startTimeEpoch);
    elapsedTimeRun = get_time_elapsed(etime(clock, startTimeRun));
    % fprintf('Time elapsed: %3.0f hours, %3.0f minutes, %3.0f seconds.\n', elapsedTimeRun(1), elapsedTimeRun(2), elapsedTimeRun(3));
    
    updatesToMake = nn.numUpdates - nn.updatesMade;
    
    % Estimate the time remaining:
    secondsPerUpdate = etime(clock, startTimeRun) / nn.updatesMade;
    secondsRemaining = updatesToMake * secondsPerUpdate;
    secondsParts = get_time_elapsed(secondsRemaining);
    % fprintf('Time remaining: %3.0f hours, %3.0f minutes, %3.0f seconds.\n', secondsParts(1), secondsParts(2), secondsParts(3));
  end
end
