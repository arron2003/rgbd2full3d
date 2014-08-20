% Returns a map of regions, where each region corresponds to a single
% object instance.
%
% Args:
%   imgObjectLabels - HxW map of pixels to object category labels where a 
%                     value of 0 indicates a missing label.
%   imgInstances - HxW map of pixels to unique object instances.
%
% Returns:
%   imgRegions - a HxW map of regions ranging from 0 to N where N is the
%                number of unique object instances.
function imgRegions = get_regions_from_labels(imgObjectLabels, imgInstances)
  [H, W] = size(imgObjectLabels);
  imgRegions = zeros(H, W);
  
  instanceMasks = get_instance_masks(imgObjectLabels, imgInstances);
  N = size(instanceMasks, 3);
  
  for ii = 1 : N
    imgRegions(instanceMasks(:,:,ii)) = ii;
  end
end