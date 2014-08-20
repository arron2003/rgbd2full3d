%% get training data for layout classifier
load('../config/splits.mat');
load ('../config/layoutLocation.mat');
opt = struct;
opt.tolErr = 0.025; % sigma of point distribution, 95% are within 2 sigma (0.05*Z)
opt.tolPara = 0.07799; % sigma cos > tolPara, we think two planes are parallel
opt.nms = 0.15; % old school nms
opt.layoutLoc = layoutLoc;
opt.gridSize = 0.03;
opt.nBadSample = 10;
trainData = struct;
trainData.feat = {};
trainData.label = {};
trainData.location = {};

allFeat = [];
allLabel = {};
allCorrect = [];

for i=trainNdxs'
  thisFeat = [];
  thisLabel = {};
  thisCorrect = [];
  thisGTLocation = [];
  fn = dir(sprintf('../processed_data/%d_*.mat', i));
  fn = fn.name;
  fprintf('%s\n', fn);
  load(sprintf('../processed_data/%s', fn));
  load(sprintf('../pixelmap/%s', fn));
  load(sprintf('../mat/%s', fn));
  data.pm = pm;
  
  % getting positive data
  for k=1:numel(model.objects)
    if ismember(model.objects{k}.model.label, {'wall', 'floor', 'ceiling'})
      for p=1:numel(model.objects{k}.model.surfaces)
        % for each surface
        plane = model.objects{k}.model.surfaces{p}.plane;
        [label, feat, location] = getLayoutFeatures(plane, data, opt);
        % get rid of bad gt (out of scope layout planes)
        s = label;
        if strcmp(s, 'hceil')
          h = max(data.Y(:)); if ((h+opt.gridSize)<location), continue; end
        elseif strcmp(s, 'hfloor')
          h = min(data.Y(:)); if ((h-opt.gridSize)>location), continue; end
        elseif strcmp(s, 'frontwall')
          h = min(data.Z(:)); if ((h-opt.gridSize)>location), continue; end
        elseif strcmp(s, 'leftwall')
          h = min(data.X(:)); if ((h-opt.gridSize)>location), continue; end
        elseif strcmp(s, 'rightwall')
          h = max(data.X(:)); if ((h+opt.gridSize)<location), continue; end
        end
    
        if ~isempty(location) % in case the plane is slanted and not manhanttan
          thisFeat = [thisFeat, feat];
          thisLabel = [thisLabel, {label}];
          thisCorrect = [thisCorrect, true];
          thisGTLocation = [thisGTLocation, location];
        end
      end
    end
  end
  
  % getting negative data
  sampleLabel = {'hceil', 'hfloor', 'frontwall', 'leftwall', 'rightwall'};
  badFeat = [];
  badLabel = [];
  badGTLocation = [];
  badCorrect = [];
  for s=sampleLabel
    
    if strcmp(s, 'hceil')
      h = max(data.Y(:)); p = [0 -1 0 0 ]; if (h<0), continue; end
    elseif strcmp(s, 'hfloor')
      h = min(data.Y(:)); p = [0 -1 0 0 ]; if (h>0), continue; end
    elseif strcmp(s, 'frontwall')
      h = min(data.Z(:)); p = [0 0 -1 0 ]; if (h>0), continue; end
    elseif strcmp(s, 'leftwall')
      h = min(data.X(:)); p = [-1 0 0 0 ]; if (h>0), continue; end
    elseif strcmp(s, 'rightwall')
      h = max(data.X(:)); p = [-1 0 0 0 ]; if (h<0), continue; end
    end
    bad = h*rand(1, opt.nBadSample);
    % nms w gt and among themselves
    vlist = thisGTLocation(ismember(thisLabel, s));
    for v = bad
      if ~sum(abs(v - vlist)<opt.nms)
        vlist = [vlist, v];
        badLabel = [badLabel, s];
        badGTLocation = [badGTLocation, v];
        badCorrect = [badCorrect, false];
        % fake plane, for feature extraction
        p(4) = v;
        [label1, feat, v1] = getLayoutFeatures(p, data, opt);
        if ~strcmp(s, label1) || v~=v1
          error('Feature extraction error');
        end
        badFeat = [badFeat, feat];
      end
    end
  end
  allFeat = [allFeat, thisFeat, badFeat];
  allLabel = [allLabel, thisLabel, badLabel];
  allCorrect = [allCorrect, thisCorrect, badCorrect];
end

save('../cache/trainLayout.mat', 'allFeat', 'allLabel', 'allCorrect');
%% train now
addpath('../util/nn');
addpath('../util/libsvm/');
addpath('../util/liblinear/');
load('../cache/trainLayout.mat');
layoutLabel = {'hceil', 'hfloor', 'frontwall', 'leftwall', 'rightwall'};
layoutModel = {};
for s=layoutLabel
  s = s{1}; % train separately for each class
  label = allCorrect( strcmp(s, allLabel) )+1;
  feat = allFeat(:, strcmp(s, allLabel) ).^.5;
  trainData = feat'; trainLabel = label';
  t = {};
  [trainData, trainMeans] = normalize_zero_mean(trainData);
  [trainData, trainStds] = normalize_unit_var(trainData);
  t.trainMeans = trainMeans;
  t.trainStds = trainStds;
  t.model = liblineartrain(trainLabel, sparse(trainData), '-s 0');
  layoutModel{end+1} = t;
end

save('../config/layoutModel.mat', 'layoutModel', 'layoutLabel');
