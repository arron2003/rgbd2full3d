% Returns an image showing candidate support surfaces. This function
% assumes that the entire room and corresponding points/normals have
% already been aligned such that the floor plane's surface normal points
% directly up.
%
% Args:
%   points3d - Nx3 point cloud where the XZ plane is the ground plane.
%   normals - Nx3 surface normals for every point.
%   H - the height of the image.
%   W - the width of the image.
%   angularThresh - the angular threshold to use, a value between 0 and pi.
%   disparityThresh - the threshold to use when checking whether points
%                     fall into the same plane.
%
% Returns:
%   supportSurfaceMap - a map whise values run from 0 to P where P is the
%                       number of planes in the images. A value of 0
%                       indicates that the point is not a supporting
%                       surface.
function supportSurfaceMap = get_support_surfaces(points3d, normals, ...
    H, W, angularThresh, disparityThresh)
  
  DEBUG = false;

  if nargin < 5
    angularThresh = 1.2;
  end
  if nargin < 6
    disparityThresh = 0.05;
  end

  points3d = swap_YZ(points3d);
  normals = swap_YZ(normals);
  
  normals = flip_normals_towards_interior(points3d, normals);
  
  [~, elevation] = cart2sph(normals(:,1), normals(:,2), normals(:,3));
  elevation = reshape(elevation, [H, W]);
  verticalMask = elevation > angularThresh;
  
  supportSurfaceMap = zeros(H, W);
  
  % Lay down a grid over the image which we'll sample.
  [X, Y] = meshgrid(1:10:W, 1:10:H);
  sampleGrid = zeros(H, W);
  sampleGrid(sub2ind([H, W], Y, X)) = 1;
  sampleGrid = sampleGrid .* verticalMask;
  
  supportSurfaceCount = 0;
  
  while nnz(sampleGrid) > 0
    if DEBUG
      figure(1);
      clf;
      
      subplot(1,3,1);
      imagesc(sampleGrid);
      
      subplot(1,3,2);
      imagesc(verticalMask);
      
      subplot(1,3,3);
      imagesc(supportSurfaceMap);
      pause(0.3);
    end

    % Randomly sample the grid space.
    remainingGridInds = find(sampleGrid);
    numGridInds = numel(remainingGridInds);
    
    % Select a random grid element (index into the HxW image).
    gridInd = remainingGridInds(randi(numGridInds));
    verticalMask(gridInd) = 0;

    % Grab the reference points normal and absolute location.
    refNormal = normals(gridInd, :);
    refPoint = points3d(gridInd, :);
    
    candidatePoints = points3d(verticalMask, :);
    numCandidates = size(candidatePoints, 1);
    
    % Take the current coordinate and project all points onto its surface
    % normal.
    disparities = (candidatePoints - repmat(refPoint, [numCandidates 1])) * refNormal';
    
    % Figure out which of these points belong to the current support
    % surface.
    vv = verticalMask;
    [Y, X] = ind2sub([H, W], find(vv));
    
    Y(abs(disparities) > disparityThresh) = [];
    X(abs(disparities) > disparityThresh) = [];
    
    iii = sub2ind([H W], Y, X);
    
    supportSurfaceMask = false(H * W, 1);
    supportSurfaceMask(iii) = 1;
    
    supportSurfaceCount = supportSurfaceCount + 1;
    supportSurfaceMap(supportSurfaceMask) = supportSurfaceCount;
    
    
    verticalMask(supportSurfaceMask) = 0;
    sampleGrid = sampleGrid & verticalMask;
  end
end
