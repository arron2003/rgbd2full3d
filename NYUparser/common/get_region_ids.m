function [regionIds, R] = get_region_ids(imgRegions)
  regionIds = unique(imgRegions);
  regionIds(regionIds == 0) = [];
  R = numel(regionIds);
end