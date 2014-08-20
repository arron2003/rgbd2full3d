% Gets a sample from the given 3D point cloud.
%
% Args:
%   points3d - the 3D point cloud to sample from, an Nx3 matrix.
%   M - the maximum number of points to return.
%   pcntInliers - the percentage of inliers from which to select from.
%
% Returns:
%   points3d - the sampled 3D point cloud.
%   inds - the indices from the original 3D point cloud that were selected.
function [points3d, inds] = get_pcd_sample(points3d, M, pcntInliers)
  if nargin < 2
    M = 100;
  else
    assert(M > 0, 'The number of points to return must be positive');
  end
  
  if nargin < 3
    pcntInliers = .9;
  else
    assert(pcntInliers <= 1, 'pcntInliers must be a percentage (0-1).');
    assert(pcntInliers >= 0, 'pcntInliers must be a percentage (0-1).');
  end
  
  assert(size(points3d, 1) >= 3);
  assert(size(points3d, 2) == 3);
  
  sampleSize = min([M, size(points3d, 1)]);
  seq = randperm2(size(points3d, 1), sampleSize);
  points3d = points3d(seq, :);
  centroid = mean(points3d);
  %dists0 = distSqr(double(centroid)', double(points3d)');
  dists = pdist2(double(centroid), double(points3d));
  [~, inds] = sort(dists, 'ascend');
  numInliers = ceil(size(points3d, 1) * pcntInliers);
  
  inds = inds(1:numInliers);
  points3d = points3d(inds, :);
end
