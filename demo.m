% RGBD to 3D demo
setup;

opt = struct('name', '1_kitchen_0004_1', 'debug', 1, 'grid_size', 0.015, 'nms', 0.15,...
  'maxNRegion', 100, 'too_small', 500, 'maxRetrieval', 10, 'missing_C', .3, 'complete_C', 3, ...
  'povray_bin', '/home/guo29/local/bin/povray' );

load config/splits.mat
addpath('preprocess/');
addpath(genpath('util/'));
addpath(genpath('NYUparser/'));
addpath(genpath('layoutDetect/'));
addpath(genpath('objectDetect/'));
addpath(genpath('mesaRender/'));
addpath(genpath('povray/'));

RandStream.setDefaultStream(RandStream.create('mrg32k3a', 'Seed', 0));
warning off;
data=parfor_load(['processed_data/' opt.name '.mat'], 'data'); 
data = preprocessRGBD(data.images, data.depths, data.rawDepths, opt);
save(['cache/' opt.name '.mat'], 'data'); 

load(['cache/' opt.name '.mat']); data=NYUparse(data, opt);
save(['cache/' opt.name '.mat'], 'data');

load(['cache/' opt.name '.mat']); 
data = findLayout(data, opt);
save(['cache/' opt.name '.mat'], 'data');

if 0 % using GT segmentation
  load(['cache/' opt.name '.mat']); 
  data1=load(['processed_data/' opt.name '.mat']); 
  data.gt=data1.data.gt; data.label_names = data1.data.label_names;
  data=findObjects(data, opt); 
  save(['cache/' opt.name '.mat'], 'data');

  data=estimateSceneGT(data, opt); 
  save(['cache/' opt.name '.mat'], 'data');

  r=renderPovrayColored(data, opt); data.povray = r;
  save(['cache/' opt.name '.mat'], 'data');
else % using generated segmentation
  load(['cache/' opt.name '.mat']); 
  data=autosegRGBD(data, opt); 
  save(['cache/' opt.name '.mat'], 'data');

  load(['cache/' opt.name '.mat']); 
  if isfield(data, 'gt'), data=rmfield(data, 'gt'); end
  data=findObjects(data, opt); 
  save(['cache/' opt.name '.mat'], 'data');

  load(['cache/' opt.name '.mat']); 
  data=estimateSceneProp(data, opt); 
  save(['cache/' opt.name '.mat'], 'data');

  r=renderPovrayColored(data, opt); data.povray = r;
  save(['cache/' opt.name '.mat'], 'data');
end

