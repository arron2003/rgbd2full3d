% Rotates the 3D points, surface normals and plane equations.
%
% Args:
%   R - 3x3 rotation matrix.
%   X - HxW matrix of X variables.
%   Y - HxW matrix of Y variables.
%   Z - HxW matrix of Z variables.
%   normals - Nx3 matrix of surface normals where N=Hx3
%   planeEqsf - Px4 matrix of plane equations.
%
% Returns:
%   X - HxW matrix of rotated X variables.
%   Y - HxW matrix of rotated Y variables.
%   Z - HxW matrix of rotated Z variables.
%   normals - Nx3 matrix of rotated surface normals where N=Hx3
%   planeEqs - Px4 matrix of plane equations.
function [X, Y, Z, normals, planeEqs] = ...
    rotate_scene(R, X, Y, Z, normals, planeEqs)

  points3d = [X(:) Y(:)  Z(:)];
  X = reshape(points3d * R(1, :)', size(X));
  Y = reshape(points3d * R(2, :)', size(X));
  Z = reshape(points3d * R(3, :)', size(X));

  normals = normals * R';

%   t = [0 ; 0 ; 0];
%   K = [K.fx 0 K.u0 ; 0 K.fy K.v0 ; 0 0 1];
%   Pf = [K*R' t]; % maps from 3d points to image

  planeEqs = [(R*planeEqs(:, 1:3)')' planeEqs(:, 4)];
end