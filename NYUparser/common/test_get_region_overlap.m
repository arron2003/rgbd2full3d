function test_suite = test_get_region_overlap
  initTestSuite;
end

function setup()
  addpath('..');
end

function test_overlap()
  imgRegions = get_regions(1);
  imgRegionsPred = imgRegions;
  R = max(imgRegions(:));
  
  N = numel(imgRegions);
  
  % Swap a random series of indices.
  inds = randi(N, [1000 2]);
  
  for ii = 1 : 1000
    tmp = imgRegionsPred(inds(ii,1));
    imgRegionsPred(inds(ii,1)) = imgRegionsPred(inds(ii,2));
    imgRegionsPred(inds(ii,2)) = tmp;
  end
  
  overlap = get_region_overlap(imgRegions, imgRegionsPred);
  for ii = 1 : R
    for jj = 1 : R
      isect = nnz(imgRegions == ii & imgRegionsPred == jj);
      union = nnz(imgRegions == ii | imgRegionsPred == jj);
      assertEqual(isect / union, overlap(ii,jj));
    end
  end
end