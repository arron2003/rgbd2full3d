% Infers the maximum extent of each plane.
%
% Args:
%   planes - Px4 matrix of P plane equations.
%   imgLabels - HxW label map whose values range from 0 (no plane) to P
%               where P is the total number of planes found.
%   points3d - Nx3 matrix of 3D world (point cloud) points.
%
% Returns:
%   planeMaps - HxWxP final plane masks.
function planeMaps = extrapolatePlanesGC(planes, imgLabels, points3d)

  % Turn on to visualize intermediate results.
  DEBUG = false;

  %%%%%%%%%%%%%%
  % parameters %
  %%%%%%%%%%%%%%
  
  % cost of assigning pixel with plane label in imgLabels to that plane
  labelUnary = -50;
  
  % multiplier for unary term based on distance
  distUnary = 50;
  
  % multiplier for unary term based on pairwise
  distPair = 100;
  
  % sets neighborhood to include diagonal edges
  neighborhood = 1;
  
  % initialize
  [H, W] = size(imgLabels);
  P = size(planes, 1);
  if size(points3d, 2) == 3
    points3d = [points3d ones(H*W, 1)];
  end
  
  planeMaps = false(H, W, P);

  for pp = 1 : P
    % Compute distance to plane and set sign so that a point has a positive 
    % value iff it is on the opposite side of the plane from the camera
    dist = points3d * planes(pp, :)';    
    if planes(pp, 4) > 0
      dist = -dist; 
    end
    dist = reshape(dist, [H W]);

    if ~any(imgLabels == pp)
      continue;
    end

    % Set unary (postive value indicates high cost to plane assignment):
    %   Points initially assigned to plane should prefer positive labels
    %   Other points with dist<0 (occluders) should have no preference
    %   Other points with dist>0 should prefer to be negative  
    unary = zeros(size(imgLabels));
    unary(imgLabels == pp) = labelUnary;
    ind = imgLabels ~= pp;
    unary(ind) = max(dist(ind) * distUnary, 0);
    
    %unary(:, [1 end]) = 1; % these would prevent flooding to borders
    %unary([1 end], :) = 1;

    % Set pairwise:
    %   Cuts between plane pixel and foreground pixel have high cost
    %   Cuts between plane pixel and background pixel have no cost
    pairwise = double(max(-dist * distPair, 0));
    pairwise = max(pairwise, 0.1);
    pairwise(:) = 1;

    unary = cat(3, unary, zeros(H, W));
    [gcseg, gc_energy] = mex_graphCut_fillPlane(-unary, pairwise, neighborhood);
    planeMaps(:,:,pp) = logical(gcseg);
    
    if DEBUG
      figure(100);
      clf;
      imagesc(planeMaps(:,:,pp));
      title(sprintf('Plane Map (Inferred %d/%d planes). Hit any key.', pp, P));
      pause;
    end
  end
end