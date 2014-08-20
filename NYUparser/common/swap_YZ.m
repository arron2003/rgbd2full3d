% Swaps the Y and Z axis in the point cloud. This function is necessary
% because the point clouds are stored in the image-coordinate frame where
% the Y axis points up (so the XY axes span the image plane) whereas many
% of the features are calculated in the right handed cartesian coordinate
% space where the Z axis points up.
%
% Args:
%   points3d - Nx3 matrix of 3D points where the Y axis points up.
%
% Returns:
%   points3d - Nx3 matrix of 3D points where the Z axis points up.
function points3d = swap_YZ(points3d)
  tmp = points3d;
  points3d(:,2) = tmp(:,3);
  points3d(:,3) = tmp(:,2);
end