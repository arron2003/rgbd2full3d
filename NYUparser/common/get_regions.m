% Loads the region map indicated by the parameters and image number.
%
% Args:
%   ii - the image number.
%   params - the parameters struct.
%
% Returns:
%   imgRegions - the region map induced by the labels if the region source
%                is consts.REGION_SRC_LABELS or the region map inferred by
%                the bottom up segmentation if the region source is
%                consts.REGION_SRC_BOTTOM_UP.
function imgRegions = get_regions(ii, params)
  Consts;

  if nargin < 2
    params = struct();
    params.regionSrc = consts.REGION_SRC_LABELS;
    params.stage = 0;
  end

  switch params.regionSrc
    case consts.REGION_SRC_LABELS
     load(sprintf(consts.imageRegionsFilename, ii), 'imgRegions');
    case consts.REGION_SRC_BOTTOM_UP
      if params.stage < 3 && ...
          (params.seg.featureSet == consts.BFT_RGBD_SUP || ...
           params.seg.featureSet == consts.BFT_RGBD_SUP_SC)
        filename = sprintf(consts.boundaryInfoPostMerge, ...
          consts.BFT_RGBD, params.stage, ii);
      else
        %filename = sprintf(consts.boundaryInfoPostMerge, ...
        %  params.seg.featureSet, params.stage, ii);
        filename = sprintf(consts.boundaryInfoPostMerge, ...
          3, params.stage, ii);
      end
      load(filename, 'boundaryInfo');
      imgRegions = boundaryInfo.imgRegions;
    otherwise
      error('not implemented...');
  end
end
