% Returns masks for the interior and border components of each region in the given boundaryInfo 
% struct.
%
% Args:
%   boundaryInfo - struct containing boundary data.
%   sz (optional) - the size (in pixels) of the boundary region.
%
% Returns:
%   regionInteriors - HxWxR set of binary masks representing the interior components of each image
%                     region.
%   regionBorders - HxWxR set of binary masks representing the border components of each image
%                     region.
%   imgBoundary - HxW image of boundary regions.
function [regionInteriors, regionBorders, imgBoundary] = ...
    get_region_components(boundaryInfo, sz)
  
  error(nargchk(1,2,nargin));
  
  if nargin < 2
    sz = 3;
  end

  se = strel('disk', sz, 8);
  
  imgBoundary = true(size(boundaryInfo.imgRegions));
  for ii = 1 : boundaryInfo.ne
    imgBoundary(boundaryInfo.edges.indices{ii}) = false;
  end
  imgBoundary2 = imerode(imgBoundary, se);
  
  [regionInteriors, regionBorders] = mex_get_region_edge_and_interior(...
      int32(boundaryInfo.imgRegions), max(boundaryInfo.imgRegions(:)), imgBoundary2);
end
