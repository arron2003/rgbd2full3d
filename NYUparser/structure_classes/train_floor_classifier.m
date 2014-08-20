% Script for training the floor classifier, called by either:
%   - run_train_floor_classifier_gt.m
%   - run_train_floor_classifier_seg.m

nn_configure_path;
nn_consts;

% Load the train test split:
load(consts.splitsPath, 'trainNdxs', 'testNdxs');

% Load the dataset.
datasetFilename = sprintf(consts.structureFeaturesDataset, ...
      params.regionSrc, params.seg.featureSet, params.stage);
load(datasetFilename, 'trainData', 'trainLabels', ...
  'testData', 'testLabels');
trainLabels = double(trainLabels == 1);
testLabels = double(testLabels == 1);

%%
fprintf('Normalizing data ...');
[trainData, trainMeans] = normalize_zero_mean(trainData);
[trainData, trainStds] = normalize_unit_var(trainData);

testData = normalize_zero_mean(testData, trainMeans);
testData = normalize_unit_var(testData, trainStds);
fprintf('DONE\n');

D = size(trainData, 2);

%%

% lambdas = [.000005, .00001, .00005 .0001 .0005 .001 .005];
lambdas = .001;

allAccuries = zeros(numel(lambdas), 2);
allMeanDiags = zeros(numel(lambdas), 2);

% lambdas = [.00005];
% ii = 1;

for ii = 1 : numel(lambdas)
  RandStream.setDefaultStream(RandStream.create('mrg32k3a', 'Seed', 1));

  nn_consts;
  nn = nn_create(D, 'region_classifier');
  nn = nn_add_layer(nn, 1, BINARY);

  nn.eta = 0.0005;
  nn.numUpdates = 5000;
  nn.lambda = lambdas(ii);
  nn.resampleStrategy = RESAMPLE_RANDOM;
  nn.resampleRatio = 1;

  nn = nn_train_sgd(nn, trainData, trainLabels);

  % Evaluate
  [accTrain, cmTrain, ranksTrain] = nn_eval(nn, trainData, trainLabels);
  [accTest, cmTest, ranksTest] = nn_eval(nn, testData, testLabels);
  fprintf('Acc Train: %f\n', accTrain);
  fprintf('Acc Test: %f\n', accTest);

  fprintf('Mean diag (Train): %f\n', mean(diag(normalize_conf_mat(cmTrain))));
  fprintf('Mean diag (Test): %f\n', mean(diag(normalize_conf_mat(cmTest))));
  
  allAccuries(ii, 1) = accTrain;
  allAccuries(ii, 2) = accTest;
  
  allMeanDiags(ii, 1) = mean(diag(normalize_conf_mat(cmTrain)));
  allMeanDiags(ii, 2) = mean(diag(normalize_conf_mat(cmTest)));
end

%% Save the results and the normalization variables to disk.
fprintf('Saving classifiers...');
outFilename = sprintf(consts.floorClassifier, params.regionSrc, params.stage);
save(outFilename, 'nn', 'cmTrain', 'cmTest', 'trainMeans', 'trainStds', '-v7.3');
fprintf('DONE.\n');