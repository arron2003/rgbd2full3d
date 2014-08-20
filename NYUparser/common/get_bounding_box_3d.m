% Returns a 3D bounding box for the given set of points. The bounding box
% is aligned with the vertical (Z) axis as is calculated by finding the
% minimum bounding rectangle in the XY (floor) plane and then finding the
% min and max values in height (Z).
%
% Args:
%   points3d - a 3D point cloud, an Nx3 matrix of points.
%   sample - (optional) whether or not to take a sample of the points for
%            calculating the bounding box.
%
% Returns:
%   boundingBox - struct which contains a basis and coefficients for the 3D
%                 bounding box of the given point cloud.
function boundingBox = get_bounding_box_3d(points3d, sample)
  if nargin < 2
    sample = true;
  else
    assert(isa(sample, 'logical'));
  end
  
  % Subsample the point cloud for two reasons:
  %   1. To get rid of outliers
  %   2. Speed
  if sample
    points3d = get_pcd_sample(points3d, 200, .9);
  end
  
  boundingBox = struct();

  % Let one basis be the up-down vector. Rows are vectors.
  boundingBox.basis = zeros(3,3);
  boundingBox.basis(3,:) = [0 0 1];

  % Project everything to the 2d plane.
  bb2d = get_min_bounding_rectangle(points3d(:,1:2), 10);
  boundingBox.basis(1:2,1:2) = bb2d.basis;
  
  maxHeight = max(points3d(:,3));
  minHeight = min(points3d(:,3));
  
  midPointZ = (maxHeight - minHeight) / 2;
  boundingBox.centroid = [bb2d.centroid midPointZ + minHeight];
  boundingBox.coeffs = [bb2d.coeffs midPointZ];
end
