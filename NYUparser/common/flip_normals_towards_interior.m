% Flips a series of surface normals so that they point away from the
% interior of the surface of the object being modeled (and towards,
% effectively, the interior of the room).
%
% This normalization ASSUMES that the origin, (0,0,0) is at the camera and
% all points are effectively IN FRONT of the viewer.
%
% Args:
%   basePoints - an Nx3 matrix of (X,Y,Z) coordinates.
%   normals - an Nx3 matrix of surface normals.
%
% Returns:
%   normals - an Nx3 matrix of surface normals, which direction has been
%             flipped if it was previously pointing 'into the surface'.
function normals = flip_normals_towards_interior(basePoints, normals)
  X = normals(:,1);
  Y = normals(:,2);
  Z = normals(:,3);

  unitBasePoints = basePoints ./ repmat(sqrt(sum(basePoints.^2,2)), [1, 3]);

  thetas = acos(sum(-unitBasePoints .* normals,2));
  flipInds = thetas > pi/2;
  
  
  % The heuristic above will fail if the point is very near the origin.
  nearOrigin = all(basePoints < 10e-8,2);
  
  % What about when the basePoints are orthogonal to the normals?
  equalInds = thetas == pi/2;
  
  % Then, just make sure theyre pointing towards the viewer.
  facingAway = Y > 0;
  X(flipInds | nearOrigin | (facingAway & equalInds)) = -X(flipInds | nearOrigin | (facingAway & equalInds));
  Y(flipInds | nearOrigin | (facingAway & equalInds)) = -Y(flipInds | nearOrigin | (facingAway & equalInds));
  Z(flipInds | nearOrigin | (facingAway & equalInds)) = -Z(flipInds | nearOrigin | (facingAway & equalInds));
  
  normals = [X Y Z];
end