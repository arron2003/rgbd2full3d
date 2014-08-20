addpath('../util/nn');
addpath('../util/libsvm/');
addpath('../util/liblinear/');
load ../cache/trainLayout.mat
sampleLabel = {'hceil', 'hfloor', 'frontwall', 'leftwall', 'rightwall'};
layoutThreshold = [];
for s=sampleLabel
  % train separately for each class
  s = s{1};
  label = allCorrect( strcmp(s, allLabel) )+1;
  feat = allFeat(:, strcmp(s, allLabel) ).^.5;
  
  % manual cv
  idx1 = randperm(numel(label)); idx1 = idx1(1:2:end); idx2 = setdiff(1:numel(label), idx1);
  trainData = feat(:,idx1)'; testData = feat(:,idx2)';
  trainLabel = label(idx1)'; testLabel = label(idx2)'; 
  [trainData, trainMeans] = normalize_zero_mean(trainData);
  [trainData, trainStds] = normalize_unit_var(trainData);
  testData = normalize_zero_mean(testData, trainMeans);
  testData = normalize_unit_var(testData, trainStds);
  %% try softmax nn
  C = numel(unique(trainLabel));
  D = size(trainData, 2);
  RandStream.setDefaultStream(RandStream.create('mrg32k3a', 'Seed', 1));
  lambdas = .0001;

  nn_consts;
  nn = nn_create(D, 'region_classifier');
  nn = nn_add_layer(nn, C, SOFTMAX);

  nn.eta = 0.0005;
  nn.numUpdates = 2000;
  nn.lambda = lambdas;

  nn = nn_train_sgd(nn, trainData, trainLabel);

  % Evaluate nn
  [accTrain, cmTrain, ranksTrain] = nn_eval(nn, trainData, trainLabel);
  [accTest, cmTest, ranksTest] = nn_eval(nn, testData, testLabel);
  outputs = nn_feed_forward(nn, testData);
  %% try liblinear
  m = liblineartrain(trainLabel,sparse(trainData), '-B 1');
  [p, acc, score1] = liblinearpredict(testLabel,sparse(testData), m);
  score1 = score1*(2*(m.Label(1)==2)-1);
  %% try libsvm
  m = svmtrain(trainLabel, trainData, '-t 4');
  [p, acc, score2] = svmpredict(testLabel,testData, m);
  score2 = score2*(2*(m.Label(1)==2)-1);
  %% plot results
  figure(1); clf;
  [X,Y,T] = perfcurve(testLabel, score1, 2, 'xCrit', 'reca', 'yCrit', 'prec');
  T(min(find(X>.5)))
  layoutThreshold = [layoutThreshold, T(min(find(X>.5)))];
  plot(X,Y,'b-');
  hold on;
  [X,Y,T] = perfcurve(testLabel, outputs(:,2), 2, 'xCrit', 'reca', 'yCrit', 'prec');
  plot(X,Y,'r-');
  [X,Y,T] = perfcurve(testLabel, score2, 2, 'xCrit', 'reca', 'yCrit', 'prec');
  plot(X,Y,'g-');
  
  
  pause;
  
end
