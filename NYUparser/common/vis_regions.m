% Visualizes the given set of regions by labeling each one with the given
% vector of values.
%
% Args:
%   imgRegions - HxW map of regions, ranging from 0 (unlabled region) to R
%                where R is the number of regions in the image.
%   values - Rx1 vector of values where each value is used to 'paint-in'
%            the corresponding region.
function img = vis_regions(imgRegions, values)
  R = max(imgRegions(:));
  
  assert(numel(values) == R);
  
  img = zeros(size(imgRegions));
  for rr = 1 : R
    img(imgRegions == rr) = values(rr);
  end
  
  if nargout == 0
    imagesc(img);
  end
end
