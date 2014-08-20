% Converts cartesian coordinates to spherical ones.
%
% Args:
%   X - X coordiantes, an Nx1 vector
%   Y - Y coordiantes, an Nx1 vector
%   Z - Z coordiantes, an Nx1 vector
% 
% Returns:
%   radius - radius of the point
%   theta - inclincation
%   psi - azimuth on the XY plane
function [radius, theta, psi] = cart2sphere(X, Y, Z)
  radius = sqrt(X.^2 + Y.^2 + Z.^2);
  theta = acos(Z ./ (radius + eps));
  psi = atan2(Y, X);
  psi(psi < 0) = 2*pi+psi(psi < 0);
end
