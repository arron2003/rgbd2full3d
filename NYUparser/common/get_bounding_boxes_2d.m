% Returns a struct containing bounding boxes for each region in the given
% region map.
%
% Args:
%   imgRegions - the region map, an HxW image where each value ranges from
%                0 (indicating a missing region) to R, the number of
%                regions in the image.
%
% Returns:
%   boundingBoxes - struct array of bounding boxes.
function boundingBoxes = get_bounding_boxes_2d(imgRegions)
  [regionIds, R] = get_region_ids(imgRegions);
  bounds = mex_bounding_boxes_2d(int32(imgRegions), int32(regionIds));
  bounds = bounds';
  
  assert(size(bounds,1) == R);
  
  boundingBoxes = struct();
  
  for ii = 1 : R
    boundingBoxes(ii).minX = bounds(ii,1);
    boundingBoxes(ii).maxX = bounds(ii,2);
    boundingBoxes(ii).minY = bounds(ii,3);
    boundingBoxes(ii).maxY = bounds(ii,4);
    
    boundingBoxes(ii).height = ...
        boundingBoxes(ii).maxY - boundingBoxes(ii).minY;
    boundingBoxes(ii).width = ...
        boundingBoxes(ii).maxX - boundingBoxes(ii).minX;
    
    boundingBoxes(ii).centroid = ...
        [boundingBoxes(ii).minX + boundingBoxes(ii).width / 2 ...
         boundingBoxes(ii).minY + boundingBoxes(ii).height / 2];
  end
end