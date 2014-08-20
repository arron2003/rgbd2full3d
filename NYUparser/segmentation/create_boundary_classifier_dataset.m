% Creates training and test sets from the boundary features computed at the current stage.
%
% Args:
%   stage - the current stage of boundary classification/merging.
%   trainNdxs - the training indicies.
%   featureSet - the parameter indicating the featureSet in use (See Consts.BFT*)
%
% Returns:
%   trainData - the training boundary feature descriptors, a NxD matrix where N is the number of
%               data points and D is the dimensionality of the boundary features.
%   testData - the testing boundary feature descriptors, a MxD matrix where M is the number of test
%              data points and D is the dimensionality of the boundary features.
%   trainLabels - Nx1 vector of labels (1,-1) indicating whether each corresponding boundary feature
%                 descriptor indicates a true boundary or not.
%   testLabels - Mx1 vector of labels (1,-1) indicating whether each corresponding boundary feature
%                descriptor indicates a true boundary or not.
function [trainData, testData, trainLabels, testLabels] = ...
    create_boundary_classifier_dataset(stage, trainNdxs, featureSet)
  Consts;
  
  switch featureSet
    case consts.BFT_RGB
      D = 13;
    case consts.BFT_D
      D = 38;
    case consts.BFT_RGBD
      D = 51;
    case consts.BFT_RGBD_SUP
      if stage < 3
        D = 51;
      else
        D = 55;
      end
    case consts.BFT_RGBD_SUP_SC
      if stage < 3
        D = 51;
      else
        D = 61;
      end
  end
  
  numTrain = numel(trainNdxs);
  numTest = consts.numImages - numTrain;

  % Initialize the training and test sets to speed up computation slightly.
  trainData = zeros(numTrain * 1000, D);
  testData = zeros(numTest * 1000, D);
  
  trainLabels = zeros(numTrain * 1000, 1);
  testLabels = zeros(numTest * 1000, 1);
  
  trainStart = 1;
  testStart = 1;
  
  fprintf('Loading all boundary features from disk:\n');
  for ii = 1 : consts.numImages
    
    fprintf('Loading boundary features (%d/%d)\r', ii, consts.numImages);
    if ~consts.useImages(ii)
      continue;
    end
    
    boundaryFeaturesFilename = sprintf(consts.boundaryFeaturesFilename, ...
        featureSet, stage, ii);
    load(boundaryFeaturesFilename, 'boundaryFeatures', 'boundaryLabels');
    
    assert(~any(isnan(boundaryFeatures(:))));
    numFeatures = numel(boundaryLabels);
    
    if isin(ii, trainNdxs)
      trainEnd = trainStart + numFeatures - 1;
      trainData(trainStart:trainEnd,:) = boundaryFeatures;
      trainLabels(trainStart:trainEnd) = boundaryLabels;
      trainStart = trainStart + numFeatures; 
    else
      testEnd = testStart + numFeatures - 1;
      testData(testStart:testEnd,:) = boundaryFeatures;
      testLabels(testStart:testEnd) = boundaryLabels(:);
      testStart = testStart + numFeatures;
    end
  end
  
  fprintf('\n');
  fprintf('================================================\n');
  fprintf('Finished loading boundary features for stage %d!\n', stage);
  fprintf('================================================\n');
  
  trainData = trainData(1:trainEnd,:);
  trainLabels = trainLabels(1:trainEnd);
  
  testData = testData(1:testEnd,:);
  testLabels = testLabels(1:testEnd);
end
