% Contains the constants used throughout the project.
consts = struct();

% Where the labeled data is stored.
consts.datasetDir = './complete/';

% Absolute path to the dataset.
consts.datasetPath = [consts.datasetDir 'nyu_depth_v2_labeled.mat'];

% The absolute path to the train/test split.
consts.splitsPath = [consts.datasetDir 'splits.mat'];

% The absolute path to the SPAMS framework:
%   http://spams-devel.gforge.inria.fr/
consts.spamsPath = './spams/spams-matlab/build/';

% Whether or not to use Gurobi. If useGurobi is false, then Matlab's
% built-in LP solver (linprog) will be used.
consts.useGurobi = true;

% The absolute path to Gurobi:
%   http://www.gurobi.com
consts.gurobiPath = '~/code1/gurobi/latest/linux64/matlab/';

% Absolute path to the Support Labels.
consts.supportLabels = [consts.datasetDir '/support_labels.mat'];

% The total number of images in the dataset.
consts.numImages = 1449;

consts.useImages = true(consts.numImages, 1);
%consts.useImages = false(size(consts.useImages));
%consts.useImages(1:200) = true;
%consts.useImages(909:1200) = true;

%%%%%%%%%%
% Images %
%%%%%%%%%%

% Directory in which the RGB images are stored.
consts.imageRgbDir = [consts.datasetDir 'images_rgb/'];

% Filename for each RGB image. To use: 
%   load(sprintf(consts.imageRgbFilename, imageNum), 'imgRgb');
consts.imageRgbFilename = [consts.imageRgbDir 'rgb_%06d.mat'];

% Directory in which the depth images are stored.
consts.imageDepthDir = [consts.datasetDir 'images_depth/'];

% Filename for each depth image. To use:
%   load(sprintf(consts.imageDepthFilename, imageNum), 'imgDepth');
consts.imageDepthFilename = [consts.imageDepthDir 'depth_%06d.mat'];

% Directory in which the raw depth images are stored.
consts.imageDepthRawDir = [consts.datasetDir 'images_depth_raw/'];

% Filename for each raw depth image. To use:
%   load(sprintf(consts.imageDepthRawFilename, imageNum), 'imgDepthRaw');
consts.imageDepthRawFilename = [consts.imageDepthRawDir 'depth_raw_%06d.mat'];

% Directory in which the object labels are stored.
consts.objectLabelsDir = [consts.datasetDir 'labels_objects/'];

% Filename for each image's object labels. To use:
%   load(sprintf(consts.objectLabelsFilename, imageNum), ...
%       'imgObjectLabels');
consts.objectLabelsFilename = [consts.objectLabelsDir 'labels_%06d.mat'];
consts.objectLabelsSegFilename = [consts.objectLabelsDir 'labels_seg_%06d.mat'];

% Directory in which the instance labels are stored.
consts.instanceLabelsDir = [consts.datasetDir 'labels_instances/'];

% Filename for each image's instance labels. To use:
%   load(sprintf(consts.instanceLabelsFilename, imageNum), ...
%       'imgInstanceLabels');
consts.instanceLabelsFilename = [consts.instanceLabelsDir 'labels_%06d.mat'];

% Directory in which the structure labels are stored.
consts.structureLabelsDir = [consts.datasetDir 'labels_structure/'];

% Filename for each image's structure labels. To use:
%   load(sprintf(consts.structureLabelsFilename, imageNum), ...
%       'imgStructureLabels');
consts.structureLabelsFilename = [consts.structureLabelsDir 'labels_%06d.mat'];

% Directory in which the ground truth regions (segmentations) are stored.
consts.imageRegionsDir = [consts.datasetDir 'regions/'];

% Filename for each image's ground truth region map. To use:
%   load(sprintf(consts.imageRegionsFilename, imageNum), 'imgRegions');
consts.imageRegionsFilename = [consts.imageRegionsDir 'regions_from_labels_%06d.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature related constants %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

consts.surfaceNormalDataDir = [consts.datasetDir 'surface_normals/'];
consts.surfaceNormalData = [consts.surfaceNormalDataDir 'surface_normals_%06d.mat'];

consts.planeDataDir = [consts.datasetDir 'planes/'];
consts.planeDataFilename = [consts.planeDataDir 'plane_data_%06d.mat'];

consts.siftDir = [consts.datasetDir 'sift/'];
consts.siftRgbFilename = [consts.siftDir 'sift_rgb_ps%d_st%d_nm%d_%06d.mat'];
consts.siftDepthFilename = [consts.siftDir 'sift_depth_ps%d_st%d_nm%d_%06d.mat'];
consts.siftDataset = [consts.siftDir 'sift_dataset_ps%d_st%d_nm%d.mat'];
consts.siftDictionary = [consts.siftDir 'sift_sc_dict_ps%d_st%d_nm%d_K%d_lambda%2.2f.mat'];

consts.structureFeaturesDir = [consts.datasetDir 'structure_class_features/'];
consts.structureFeaturesFilename = [consts.structureFeaturesDir 'features_src%d_set%d_stage%d_%06d.mat'];
consts.structureFeaturesDataset = [consts.structureFeaturesDir 'dataset_src%d_set%d_stage%d.mat'];
consts.structureClassifier = [consts.structureFeaturesDir 'classifier_src%d_set%d_stage%d.mat'];

consts.floorClassifier = [consts.structureFeaturesDir 'floor_classifier_src%d_stg%d.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%
% Support-related files %
%%%%%%%%%%%%%%%%%%%%%%%%%

consts.supportFeaturesDir = [consts.datasetDir 'support/'];
consts.supportFeaturesGeometry = [consts.supportFeaturesDir 'features_geo_src%d_set%d_stg%d_%06d.mat'];
consts.supportFeaturesContainment = [consts.supportFeaturesDir 'features_con_src%d_set%d_stg%d_%06d.mat'];
consts.supportFeaturesHorz = [consts.supportFeaturesDir 'features_horz_src%d_set%d_stg%d_%06d.mat'];
consts.supportDataset = [consts.supportFeaturesDir 'dataset_src%d_set%d_stg%d.mat'];
consts.supportClassifier = [consts.supportFeaturesDir 'classifier_src%d_set%d_stg%d_res%d_rat%d.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmentation Constants %
%%%%%%%%%%%%%%%%%%%%%%%%%%

consts.watershedDir = [consts.datasetDir 'watershed/'];
consts.watershedFilename = [consts.watershedDir 'watershed_%06d.mat'];

consts.boundaryFeaturesDir = [consts.datasetDir 'boundary_features/'];
consts.boundaryFeaturesFilename = [consts.boundaryFeaturesDir 'type%d_stg%d_%06d.mat'];
consts.boundaryFeaturesDataset = [consts.boundaryFeaturesDir 'dataset_type%d_stg%d.mat'];
consts.boundaryInfoPostMerge = [consts.boundaryFeaturesDir 'info_type%d_stg%d_%06d.mat'];
consts.boundaryClassifierFilename = [consts.boundaryFeaturesDir 'classifier_type%d_stg%d.mat'];
consts.boundaryClassifierPlotFilename = [consts.boundaryFeaturesDir 'classifier_roc_type%d_stg%d.png'];

%%%%%%%%%%%
% Results %
%%%%%%%%%%%

% Directory in which support inference results are stored.
consts.resultsDir = [consts.datasetDir 'results/'];

% Filename for results of support inference using Baseline #1: Image Plane
% Rules. To use:
%   load(sprintf(consts.resultsImgFilename, imageNum), ...
%       'supportLabelsPred');
consts.resultsImgFilename = [consts.resultsDir 'src%d_img_%06d.mat'];

% Filename for results of support inference using Baseline #2: Structure
% Class Rules. To use:
%   load(sprintf(consts.resultsStrFilename, imageNum), ...
%       'supportLabelsPred', 'M');
consts.resultsStrFilename = [consts.resultsDir 'src%d_str_%06d.mat'];

% Filename for results of support inference using Baseline #3: Support
% Classifier. To use:
%   load(sprintf(consts.resultsSupFilename, imageNum), ...
%       'supportLabelsPred');
consts.resultsSupFilename = [consts.resultsDir 'src%d_sup_%06d.mat'];
consts.resultsLpFilename = [consts.resultsDir 'src%d_lp_%06d.mat'];
consts.resultsIpFilename = [consts.resultsDir 'src%d_ip_%06d.mat'];

% Regions may be obtained via the ground truth object and instance labels
% or from a bottom up segmentation.
consts.REGION_SRC_LABELS = 0;
consts.REGION_SRC_BOTTOM_UP = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary Feature Type constants are used to select different sets of
% features when performing segmentation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use RGB-features only.
consts.BFT_RGB = 1;

% Use depth-features only.
consts.BFT_D = 2;

% Use both RGB and depth features.
consts.BFT_RGBD = 3;

% Use RGB, depth and support features.
consts.BFT_RGBD_SUP = 4;

% Use RGB, depth, support features and structure class features.
consts.BFT_RGBD_SUP_SC = 5;

consts.FLOOR = 1;
consts.STRUCTURE = 2;
consts.FURNITURE = 3;
consts.PROP = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Support Inference Types %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
consts.SUP_INF_IMG_PLN_RULES = 1;
consts.SUP_INF_STR_CLS_RULES = 2;
consts.SUP_INF_LCL_CLASSIFIER = 3;
consts.SUP_INF_LP = 4;
consts.SUP_INF_IP = 5;


% If a region has fewer than this number of pixels, just drop it. This is used primarily in creating
% the final 'ground truth' regions. See run_extract_regions_from_labels.m.
consts.MIN_PIXELS_PER_REGION = 10;

