% Returns the object/class labels for each of the segments in the proposed
% region map.
%
% Args:
%   imgRegions - HxW image of regions 1..R
%   instanceMasks - HxWxM where M is the number of instances and each layer
%                   contains a binary mask for each instance.
%   instanceClasses - the class identities for each of the masks.
%
% Returns:
%   segClassLabels - class labels for each proposed segment.
%   segInstanceLabels - instance labels for each poposed segment.
%   purity - the percentage of the true mask occupied by each segment.
%   completeness - the percentage of the true mask made up by each segment.
function [segClassLabels, segInstanceLabels, pcntOfInstances, purity, completeness] = ...
    regions2labels(imgRegions, instanceMasks, instanceClasses)

  [~, R] = get_region_ids(imgRegions);

  missing = ~any(instanceMasks, 3);
  instanceMasks = cat(3, instanceMasks, missing);
  instanceClasses(end+1) = 0;

  numInstances = size(instanceMasks, 3);
  
  % Captures the proportion each segment occupies out of each instance map.
  pcntOfInstances = zeros(R, numInstances); 
  segInstanceLabels = zeros(R, 1);
  segClassLabels = zeros(R, 1);
  purity = zeros(R, 1);
  completeness = zeros(R, 1);

  stats = regionprops(imgRegions, 'PixelIdxList', 'Area');
  indsPerSegment = {stats.PixelIdxList};
  numPixelsPerSegment = [stats.Area];

  pixelsPerInstance = reshape(sum(sum(instanceMasks, 1), 2), [1 numInstances]);

  for ii = 1 : numInstances
    instanceMask = instanceMasks(:,:,ii);
    for rr = 1 : R
      pcntOfInstances(rr, ii) = nnz(instanceMask(indsPerSegment{rr})) / numPixelsPerSegment(rr);
    end
  end

  assert(all(abs(sum(pcntOfInstances,2) - 1)< 10e-5));
  
  for rr = 1 : R

    % If at least 50% of the pixels overlap with two separate instancs
    if sum((segInstanceLabels(rr, :) > 0.5) >= 2) 
      % Choose one with better intersection over union.
      [~, mi] = max(pcntOfInstances(rr, :) * numPixelsPerSegment(rr) ./ pixelsPerInstance);
    else
      % Choose the majority label.
      [~, mi] = max(pcntOfInstances(rr, :));
    end
    
    segClassLabels(rr) = instanceClasses(mi);
    segInstanceLabels(rr) = mi;
    
    % If the label is missing...
    if segClassLabels(rr) == 0
      purity(rr) = pcntOfInstances(rr, end);
    else
      % ignore the void labels
      purity(rr) = max(pcntOfInstances(rr, 1:end-1)) / (1-pcntOfInstances(rr, end));
    end
    
    completeness(rr) = pcntOfInstances(rr, mi) * numPixelsPerSegment(rr) / pixelsPerInstance(mi);
  end
end
