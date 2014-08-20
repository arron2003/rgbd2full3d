% Projects the depth points from the image plane to the 3D world
% coordinates.
%
% Args:
%   imgDepth - depth map which has already been projected onto the RGB
%              image plane, an HxW matrix where H and W are the height and
%              width of the matrix, respectively.
%
% Returns:
%   points3d - the point cloud in the world coordinate frame, an Nx3
%              matrix.
function points3d = rgb_plane2rgb_world(imgDepth)
  camera_params;
  [H, W] = size(imgDepth);

  % Make the original consistent with the camera location:
  [xx, yy] = meshgrid(1:W, 1:H);

  x3 = (xx - cx_rgb) .* imgDepth / fx_rgb;
  y3 = (yy - cy_rgb) .* imgDepth / fy_rgb;
  z3 = imgDepth;
  
  points3d = [x3(:) -y3(:) z3(:)];
end