% Trains a boundary classifier at the current stage.
%
% Args:
%   stage - the current stage of boundary classification.
%   trainData - the training features, an NxD matrix.
%   trainLabels - the labels for the training set, an Nx1 vector (1,-1).
%   testData - the testing features, an MxD matrix.
%   testLabels - the labels for the testing set, an Mx1 vector (1,-1).
%   params - the parameter struct (See Params.m);
%
% Returns:
%   classifier - a trained classifier for merging segments at the current stage.
function classifier = train_boundary_classifier_dt(stage, ...
    trainData, trainLabels, testData, testLabels, params)
  assert(isscalar(stage));
  assert(isstruct(params));
  
  % Make sure the training and test sets have the same dimensionality.
  assert(size(trainData,2) == size(testData,2));
  
  Consts;

  numTraining = size(trainData,1);
  
  if params.seg.maxTrainingSize < numTraining
    sampleNdxs = randperm(numTraining);
    sampleNdxs = sampleNdxs(1:params.seg.maxTrainingSize);
    trainDataSample = trainData(sampleNdxs, :);
    trainLabelsSample = trainLabels(sampleNdxs, :);
  else
    trainDataSample = trainData;
    trainLabelsSample = trainLabels;
  end
  
  numIterations = params.seg.training.numIters(stage);
  numNodes = params.seg.training.numNodes(stage);

  stopval = 0; 
  w = []; % uniform weighting
  
  fprintf('Training boundary classifier (stage %d)\n', stage);
  classifier = train_boosted_dt_2c(trainDataSample, [], trainLabelsSample, numIterations, numNodes, stopval, w);    
  
  trainScores = test_boosted_dt_mc(classifier, trainData);
  testScores = test_boosted_dt_mc(classifier, testData);

  figure;
  trainAuc = roc(trainScores, trainLabels(:) == 1, 'b');
  hold on;
  testAuc = roc(testScores, testLabels(:) == 1, 'r');
  legend({'training', 'testing'});
  hold off;
  
  title(sprintf('Train AUC: %f, Test AUC: %f', trainAuc, testAuc));
  
  fprintf('TrainAuc: %2.1f\n', 100 * trainAuc);
  fprintf('TestAuc: %2.1f\n', 100 * testAuc);
  
  plotFilename = sprintf(consts.boundaryClassifierPlotFilename, ...
    params.seg.featureSet, stage);
  
  print('-f1', '-dpng', plotFilename);
  
  fprintf('\n');
  fprintf('===================================================\n');
  fprintf('Finished training boundary classifier for stage %d!\n', stage);
  fprintf('===================================================\n');
end


