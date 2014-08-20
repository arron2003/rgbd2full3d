% Creates a grid of sample points.
%
% Args:
%   H - height of the image.
%   W - width of the image.
%   margin - the margin to use from the sides of the image to lay down the
%            grid.
%   stride - the number of pixels between grid points.
%
% Returns:
%   sampleMask - a logical HxW mask of points to sample from.
%   GH - the grid height.
%   GW - the grid width.
function [sampleMask, GH, GW] = get_sample_grid(H, W, margin, stride)
  X = margin + 1 : stride : W - margin - 1;
  Y = margin + 1 : stride : H - margin - 1;
  [X, Y] = meshgrid(X, Y);
  gridInds = sub2ind([H, W], round(Y), round(X));
  
  sampleMask = false(H, W);
  sampleMask(gridInds) = 1;
  
  % Now re-create the grid (used later when we setup the dataset).
  [Y, X] = ind2sub([H, W], gridInds);
  GH = numel(unique(Y));
  GW = numel(unique(X));
end
