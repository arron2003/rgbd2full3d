% Merges the regions based on the output of the boundary classifier.
%
% Args:
%   boundaryInfo - the boundaryInfo struct containing the region data.
%   boundaryFeatures - PxD matrix of features where P is the number of pairs of regions being
%                      considered and D is the dimensionality of the current feature set.
%   classifier - the trained classifier (See train_boundary_classifier_dt.m)
%   stage - the current stage of region merging.
%   params - the parameters struct (See Params.m)
%
% Returns:
%   result - the resulting boundary info struct.
%   pbnd - the probability of boundaries.
function [result, pbnd] = merge_regions(boundaryInfo, boundaryFeatures, ...
      classifier, stage, params)
  Consts;
  
  scores = test_boosted_dt_mc(classifier, boundaryFeatures);
  pbnd = 1./(1+exp(-scores));

  confs = pbnd;

  probOfBoundary = 1-confs;

  maxProb = params.seg.minProbs(stage);
  minRegions = 0;

  
  
  numEdges = boundaryInfo.ne;
  R = boundaryInfo.nseg;

  % spLR(i, [1 2]) gives the left and right, resp., 
  spLR = boundaryInfo.edges.spLR;

  % Keeps track of the edge IDs between regions.
  edgeletChain = num2cell(1:numEdges);

  % regionSuperpixels is a cell array, each entry of which reprents a
  % region. Each region is built from a collection of superpixels. So if
  % superpixels 3 and 4 are merged to form region 106, then:
  %   regionSuperpixels{106} = [3, 4];
  % 
  regionSuperpixels = num2cell(1:R);
  isRegionValid = true(R, 1);
  e2chain = (1:numEdges);
  needsUpdating = false(numEdges, 1);

  % unite borders between two superpixels that are split
  reme = [];
  spLRm = zeros([R R], 'uint16');
  for ii = 1:numEdges

    % Why use this min/max? instead of just picking left and right?
    s1 = min(spLR(ii, :));
    s2 = max(spLR(ii,:));

    if spLRm(s1,s2) == 0
      spLRm(s1,s2) = ii;
    else
      % Why do these 'duplicate' edges exist?
      reme(end+1) = ii;
      k2 = spLRm(s1, s2);
      edgeletChain{k2} = [edgeletChain{k2} edgeletChain{ii}];        
      e2chain(edgeletChain{k2}) = k2;
    end
  end

  spLR(reme, :) = [];
  edgeletChain(reme) = [];
  e2chainshift = ones(numel(e2chain), 1);
  e2chainshift(reme) = 0;
  e2chainshift = cumsum(e2chainshift);
  e2chain = e2chainshift(e2chain);
  clear spLRm

  probOfBoundary = probOfBoundary(1:numEdges, 1); % probability of edge being on

  %% Iteratively merge regions

  edgeCost = zeros(numel(edgeletChain), 1);
  for ii = 1:numel(edgeletChain)
    edgeCost(ii) = min(probOfBoundary(edgeletChain{ii}));
  end

  numRemainingRegions = R;
  numRemainingEdges = numel(edgeletChain);

  iter = 0;
  while numRemainingRegions > minRegions
    iter = iter + 1;                           

%     % Find the best edge:
%     [minCost, minInd] = max(edgeCost);
%     if isempty(minInd) || minCost < maxProb
%       break;
%     end
    
    % Find a random edge.
    eligibleEdges = find(edgeCost >= maxProb);
    seq = randperm(numel(eligibleEdges));
    
    if isempty(seq)
      break;
    end

    minInd = eligibleEdges(seq(1));
 
    keep = true(numRemainingEdges, 1);
    reme = minInd;
    keep(minInd) = false;

    r1 = spLR(minInd, 1);
    r2 = spLR(minInd, 2);

    newr = r1;

    % get neighbors to left and right of r1 and r2
    ind1L = find(spLR(:, 1)==r1);    
    ind2L = find(spLR(:, 1)==r2);  
    left1 = spLR(ind1L, 2);
    left2 = spLR(ind2L, 2);    
    ind1R = find(spLR(:, 2)==r1);
    ind2R = find(spLR(:, 2)==r2);
    right1 = spLR(ind1R, 1);
    right2 = spLR(ind2R, 1);     

    spLR([ind1L ; ind2L], 1) = newr;    
    spLR([ind1R ; ind2R], 2) = newr;        

    % unite any split borders that may arise
    indL = [ind1L ; ind2L ; ind1R ; ind2R]';  
    s1 = [left1 ; left2 ; right1 ; right2]';      
    for k1 = 1:numel(s1)
      for k2 = k1+1:numel(s1)
        if s1(k1)==s1(k2) && ~any(reme==indL(k2))                
          keep(indL(k2)) = false;                                
          i1 = indL(k1);
          edgeletChain{i1} = [edgeletChain{i1} edgeletChain{indL(k2)}];
          e2chain(edgeletChain{i1}) = i1;
          
          % This is where the randomness comes in - by using the min cost,
          % rather than some other metric, the order of merging will affect
          % the ultimate result.
          edgeCost(i1) = min(probOfBoundary(edgeletChain{i1}));
          %needsUpdating(i1) = true;                                
        end
      end
    end
    needsUpdating(indL) = true;

    % remove extra edges, regions
    isRegionValid(r2) = false;        

    regionSuperpixels{newr} = [regionSuperpixels{r1} regionSuperpixels{r2}];        

    edgeletChain = edgeletChain(keep);
    spLR = spLR(keep, :);
    edgeCost = edgeCost(keep);
    needsUpdating = needsUpdating(keep);

    e2chainshift = max(cumsum(keep),1);
    e2chain = e2chainshift(e2chain);

    numRemainingRegions = numRemainingRegions - 1;
    numRemainingEdges = numel(edgeCost);
  end

  rshift = cumsum(isRegionValid);
  spLR = rshift(spLR);

  result.edgeIndices = boundaryInfo.edges.indices;
  result.edgeletChain = edgeletChain;
  result.edgeCost = edgeCost;
  result.chainLR = spLR;
  result.regions = regionSuperpixels(isRegionValid);
end

