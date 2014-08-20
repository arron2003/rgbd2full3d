% Contains the parameters to the different modules of the
% segmentation/support inference pipeline.

params = struct();

% Whether or not to display intermediate results. Note that this will
% produce many figures and additional output.
params.debug = false;

% Whether or not to overwrite existing files.
params.overwrite = 1;

% The source of regions (ground truth or bottom up segmented regions).
params.regionSrc = consts.REGION_SRC_LABELS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plane Fitting Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.planes.blockWidths = [1 3 6 9];
params.planes.relDepthThresh = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIFT Feature Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.sift.stride = 10;
params.sift.patchSize = 40;
params.sift.gridMargin = 15;
params.sift.normMethod = 1;

% Parameters for SIFT Sparse Coding.
params.sc.K = 1000;
params.sc.lambda = .1;

% Parameters for region features.
params.numLevels = 2;
params.numInclinationBins = 5;
params.numAngleBins = 8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmentation Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The number of stages to use for segmentation.
params.seg.numStages = 5;

% The maximum number of data points to use when training the boundary
% classifier.
params.seg.maxTrainingSize = 500000;

% The number of iterations to run when training the decision tree.
params.seg.training.numIters = [30 30 30 30 30];

% The maximum number of nodes to use in the decision tree.
params.seg.training.numNodes = [20 20 20 10 10];

% The minimum probabilities required for merging two regions.
params.seg.minProbs = [0.9 0.8 0.7 0.7 0.7];
  
if params.regionSrc == consts.REGION_SRC_LABELS
  params.stage = 0;
  params.seg.featureSet = 0;
else
  params.stage = 5;

  % Whether or not to use an extended set of features.
  params.seg.featureSet = consts.BFT_RGBD;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supper Inference Params %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.support.classifier.eta = 0.0001;
params.support.classifier.numUpdates = 40000;
params.support.classifier.resampleStrategy = 0;
params.support.classifier.resampleRatio = 0;
params.support.verifyEnergy = false;

params.support.coeffsGt.alpha = 2; % Data cost for S.
params.support.coeffsGt.beta = 9; % Data cost for M.
params.support.coeffsGt.kappa = 5; % Global ground consistency
params.support.coeffsGt.lambda = 2; % Structure class transition Costs.
params.support.coeffsGt.tau = 1; % Vertical Consistency.
params.support.coeffsGt.upsilon = 2; % Horizontal Consistency.
params.support.coeffsGt.rho2 = 65;

% params.support.coeffsSeg.alpha = 17; % Data cost for S.
% params.support.coeffsSeg.beta = 10; % Data cost for M.
% params.support.coeffsSeg.kappa = 1; % Global ground consistency
% params.support.coeffsSeg.lambda = 1; % Structure class transition Costs.
% params.support.coeffsSeg.tau = 1; % Vertical Consistency.
% params.support.coeffsSeg.upsilon = 8; % Horizontal Consistency.
% params.support.coeffsSeg.rho2 = 518;

params.support.coeffsSeg.alpha = 2; % Data cost for S.
params.support.coeffsSeg.beta = 9; % Data cost for M.
params.support.coeffsSeg.kappa = 5; % Global ground consistency
params.support.coeffsSeg.lambda = 2; % Structure class transition Costs.
params.support.coeffsSeg.tau = 1; % Vertical Consistency.
params.support.coeffsSeg.upsilon = 2; % Horizontal Consistency.
params.support.coeffsSeg.rho2 = 500;


params.support.structPerc.numIter = 50;
