%% load stuff
addpath('../util/nn');
addpath('../util/libsvm/');
addpath('../util/liblinear/');
addpath('../annotation/');
%fn = '1_kitchen_0004_1.mat';
fn = '15_office_0008_1.mat';
load (['../processed_data/' fn]);
load (['../mat/' fn]);
load (['../pixelmap/' fn]);
load ('../config/layoutLocation.mat');
load ('../config/layoutModel.mat');
data.pm = pm;
opt = struct;
opt.tolErr = 0.025; % sigma of point distribution, 95% are within 2 sigma (0.05*Z)
opt.tolPara = 0.07799; % sigma cos > tolPara, we think two planes are parallel
opt.layoutLoc = layoutLoc;
%opt.nBadSample = 100;
opt.gridSize = 0.03;
opt.nms = 0.15; % old school nms
sampleLabel = {'hceil', 'hfloor', 'frontwall', 'leftwall', 'rightwall'};

%% start the work
model.data.coord = double(cat(3, data.X, data.Y, data.Z));
model.data.image = data.images;
clf;
opt1=get_opt(); opt1.render_pointcloud = 1; opt1.render_scale=20;
scenemodel_render(model, opt1);
  
layoutThreshold =  [0.0571    1.3018    0.4777    0.2355    0.1029]; % from .5 recall
for s=layoutLabel
  s = s{1};
  badFeat = [];
  badLabel = [];
  badLocation = [];
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
  %bad = h*rand(1, opt.nBadSample);
  bad = sign(h)*(1e-3:opt.gridSize:(abs(h)+opt.gridSize));
  for v = bad
    badLabel = [badLabel, s];
    badLocation = [badLocation, v];
    % fake plane, for feature extraction
    p(4) = v;
    [label1, feat, v1] = getLayoutFeatures(p, data, opt);
    if ~strcmp(s, label1) || v~=v1
      error('Feature extraction error');
    end
    badFeat = [badFeat, feat.^.5];
  end
  % do nms
  id = find(ismember(layoutLabel, s)); lmodel = layoutModel{id};
  testData = badFeat';
  testData = normalize_zero_mean(testData, lmodel.trainMeans);
  testData = normalize_unit_var(testData, lmodel.trainStds);
  %[p, acc, score] = svmpredict(ones(size(testData,1), 1), testData, lmodel.model);
  [p, acc, score] = liblinearpredict(ones(size(testData,1), 1), sparse(testData), lmodel.model);
  score = score*(2*(lmodel.model.Label(1)==2)-1);
  nmsGrid = opt.nms/opt.gridSize;
  isLocalMax = (score == ordfilt2(score, nmsGrid*2+1, ones(nmsGrid*2+1, 1)));
  isLocalMax = isLocalMax & score>layoutThreshold(id); % or change the threshold here
  
  v = bad(isLocalMax);
  % render proposed surfaces
  for h = v
    if strcmp(s, 'hceil') || strcmp(s, 'hfloor')
      p = [min(data.X(:)), min(data.X(:)) max(data.X(:)) max(data.X(:)) ;
           h, h, h, h;
           min(data.Z(:)), max(data.Z(:)) min(data.Z(:)) max(data.Z(:))];
      c = cat(3, zeros(2), ones(2), zeros(2));
    elseif strcmp(s, 'frontwall')
      p = [min(data.X(:)), min(data.X(:)) max(data.X(:)) max(data.X(:)); 
           min(data.Y(:)), max(data.Y(:)) min(data.Z(:)) max(data.Y(:));
           h, h, h, h];
      c = cat(3, zeros(2), zeros(2), ones(2));
    elseif strcmp(s, 'leftwall') || strcmp(s, 'rightwall')
      p = [h, h, h, h;
           min(data.Y(:)), min(data.Y(:)) max(data.Y(:)) max(data.Y(:)); 
           min(data.Z(:)), max(data.Z(:)) min(data.Z(:)) max(data.Z(:)) ];
      c = cat(3, ones(2), zeros(2), zeros(2));
    end
    p = reshape((model.camera.R*p)', [2 2 3])*opt1.render_scale;
    thandle = surf(p(:,:,1), p(:,:,2), p(:,:,3), c);
    alpha(thandle, 0.2);
  end
end
%print('-dpng','-r300','example.png');
keyboard;