% Returns the region overlap between the ground truth regions and the
% predicted regions.
%
% Args:
%   imgRegionsTrue - HxW map of pixels to regions with R1 total regions.
%   imgRegionsPred - HxW map of pixels to regions with R2 total regions.
%
% Returns:
%   overlap - R1xR2 matrix of overlaps.
%   intersection - R1xR2 matrix of intersection.
%   union - R1xR2 matrix of union.
function [overlap, intersection, union] = get_region_overlap(imgRegionsTrue, imgRegionsPred)
  assert(all(size(imgRegionsTrue) == size(imgRegionsPred)));
  numTrue = max(imgRegionsTrue(:));
  numPred = max(imgRegionsPred(:));
  [overlap, intersection, union] = mex_overlap(int32(imgRegionsTrue), ...
      int32(imgRegionsPred), int32(numTrue), int32(numPred));
end
