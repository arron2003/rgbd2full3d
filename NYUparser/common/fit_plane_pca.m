% Gets the surface normal for the given 3D points.
%
% Args:
%   points3d - Nx3 set of 3D points [X, Y, Z]
%
% Returns:
%   normal - 1x3 vector
%   basis - 2x3 matrix, the two basis vectors.
function [normal, basis] = fit_plane_pca(points3d)
  coeff = princomp(points3d);
  normal = coeff(:,3);
  normal = normal ./ sqrt(sum(normal.^2));
  normal = normal';
  basis = coeff(:,1:2)';
end
