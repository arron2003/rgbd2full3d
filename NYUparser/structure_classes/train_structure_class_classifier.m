% Trains the structure class classifier and saves it to disk.
%
% Args:
%   params - the parameter struct to guide training.
function train_structure_class_classifier(params)
  Consts;

  % Load the train test split.
  load(consts.splitsPath, 'trainNdxs');

  %% Load the train/test sets.
  fprintf('Loading dataset...');
  params.seg.featureSet = 0;
  datasetFilename = sprintf(consts.structureFeaturesDataset, ...
      params.regionSrc, params.seg.featureSet, params.stage);
  load(datasetFilename, 'trainData', 'trainLabels', ...
    'testData', 'testLabels');
  fprintf('DONE.\n');

  C = numel(unique(trainLabels));

  %%
  fprintf('Normalizing data ...');
  [trainData, trainMeans] = normalize_zero_mean(trainData);
  [trainData, trainStds] = normalize_unit_var(trainData);

  testData = normalize_zero_mean(testData, trainMeans);
  testData = normalize_unit_var(testData, trainStds);
  fprintf('DONE\n');

  D = size(trainData, 2);

  %%
  RandStream.setDefaultStream(RandStream.create('mrg32k3a', 'Seed', 1));
  lambdas = .0001;

  nn_consts;
  nn = nn_create(D, 'region_classifier');
  nn = nn_add_layer(nn, C, SOFTMAX);

  nn.eta = 0.0005;
  nn.numUpdates = 5000;
  nn.lambda = lambdas;

  nn = nn_train_sgd(nn, trainData, trainLabels);

  % Evaluate
  [accTrain, cmTrain, ranksTrain] = nn_eval(nn, trainData, trainLabels);
  [accTest, cmTest, ranksTest] = nn_eval(nn, testData, testLabels);
  fprintf('Acc Train: %f\n', accTrain);
  fprintf('Acc Test: %f\n', accTest);

  fprintf('Mean diag (Train): %f\n', mean(diag(normalize_conf_mat(cmTrain))));
  fprintf('Mean diag (Test): %f\n', mean(diag(normalize_conf_mat(cmTest))));


  %% Save the results and the normalization variables to disk.
  outFilename = sprintf(consts.structureClassifier, params.regionSrc, params.seg.featureSet, params.stage);
  fprintf('Saving classifier to file %s...', outFilename);
  save(outFilename, 'nn', 'cmTrain', 'cmTest', 'trainMeans', 'trainStds', '-v7.3');
  fprintf('DONE.\n');


end
