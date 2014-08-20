% Returns the planar disparity of the given set of 3D points from the plane
% with the given center and surface normal.
%
% Args:
%   points3d - Nx3 matrix of 3D points.
%   normal - 1x3 surface normal.
%   centroid - 1x3 vector, the 3D coordinate of the center of the plane.
%
% Returns:
%   sse - the sum of squared disparities over all of the points.
%   disparities - Nx1 vector of disparities of each point from the plane.
function [sse, disparities] = get_planar_disparity(points3d, normal, centroid)
  [N, D] = size(points3d);
  
  if nargin >= 3
    points3d = points3d - repmat(centroid, [N, 1]);
  end
  
  assert(D == 3);
  disparities = sum(points3d .* repmat(normal, [N, 1]), 2);
  sse = sum(disparities.^2);
end
