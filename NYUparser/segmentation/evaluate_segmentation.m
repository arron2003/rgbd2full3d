% Returns evaluation metrics for the given segmentation. The scoring is
% calculated by diving the intersection of each region and segment by their
% union.
%
% Args:
%   imgRegionsTrue - the ground truth segmentation map with R regions.
%   imgRegionsPred - the prediction segmentation map with S regions.
%
% Returns:
%   weightedScore - the average score, weighted by object size.
%   unweightedScore - the average score, unweighted.
%   regionScores - the overlap scores per TRUE region, an Rx1 vector.
%   bestSegments - indices of the best segments (out of S) that match each
%                  of the true regions, an Rx1 vector.
%   overlap - the raw overlap scores, an RxS matrix.
function [weightedScore, unweightedScore, regionScores, bestSegments, overlap] = ...
    evaluate_segmentation(imgRegionsTrue, imgRegionsPred)
  
  [trueRegionIds, numTrueRegions] = get_region_ids(imgRegionsTrue);
  
  overlap = get_region_overlap(imgRegionsTrue, imgRegionsPred);
  
  [regionScores, bestSegments] = max(overlap, [], 2);
  
  % First, calculate the area score.
  weightedScoreNum = 0;
  weightedScoreDenom = 0;
  for ii = 1 : numTrueRegions
    pixelAreaTrueRegion = nnz(imgRegionsTrue == trueRegionIds(ii));
    weightedScoreNum = weightedScoreNum + pixelAreaTrueRegion * max(overlap(ii,:));
    weightedScoreDenom = weightedScoreDenom + pixelAreaTrueRegion;
  end
  
  weightedScore = weightedScoreNum / weightedScoreDenom;
  
  unweightedScore = mean(regionScores);
  
  assert(~isnan(weightedScore));
end
