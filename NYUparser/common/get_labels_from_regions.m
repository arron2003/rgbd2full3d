% Gets the vector of labels assigned to each region in imgRegions.
%
% Args:
%   imgRegions - HxW image whose values range from 0 (ignored region) to R
%                (the maximum region value).
%   imgLabels - HxW image whose values range from 0 (missing label) to C
%               (the maximum label).
%
% Returns:
%   labels - Rx1 vector of labels, one for each of the R regions of the
%            region map.
function labels = get_labels_from_regions(imgRegions, imgLabels)
  error(nargchk(2,2,nargin));
  overlap = get_overlap_matrix(imgRegions, imgLabels+1);
  [~, labels] = max(overlap, [], 2);
  labels = labels - 1;
end

function overlap = get_overlap_matrix(imgRegionsTrue, imgRegionsPred)
  numTrue = max(imgRegionsTrue(:));
  numPred = max(imgRegionsPred(:));
  overlap = mex_overlap_matrix(int32(imgRegionsTrue), ...
    int32(imgRegionsPred), numTrue, numPred);
end