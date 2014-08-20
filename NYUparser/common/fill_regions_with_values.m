function imgNew = fill_regions_with_values(imgRegions, regionValues, regionIds)
  if nargin < 3
    regionIds = get_region_ids(imgRegions);
  end

  assert(numel(regionIds) == numel(regionValues));

  imgNew = zeros(size(imgRegions));
  for rr = 1 : numel(regionIds)
    imgNew(imgRegions == regionIds(rr)) = regionValues(rr);
  end
end
