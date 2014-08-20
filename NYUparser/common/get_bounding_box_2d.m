% Returns the bounding box of the given binary mask.
% 
% Args:
%   mask - HxW logical image where pixels that are turned on indicate the
%          shape for which we want a bounding box.
%
% Returns:
%   bb - a struct containing fields indicating the bounding boxes minimum
%        and maximum rows and columns, the height and width of the bounding
%        box, and the coordinate of its centroid.
function bb = get_bounding_box_2d(mask)
  assert(isa(mask, 'logical'));

  bb = struct();
  bb.minX = find(sum(mask, 1), 1, 'first');
  bb.maxX = find(sum(mask, 1), 1, 'last');
  bb.minY = find(sum(mask, 2), 1, 'first');
  bb.maxY = find(sum(mask, 2), 1, 'last');
  
  bb.height = bb.maxY - bb.minY;
  bb.width = bb.maxX - bb.minX;
  
  bb.centroid = [bb.minY + bb.height/ 2 ...
                 bb.minX + bb.width / 2];
end
