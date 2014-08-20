function data = autosegRGBD( data, opt )
  addpath(genpath([fileparts( mfilename('fullpath') ) '/proposals']));
  pm=data.pm(:,:,4); nSps = data.baseSeg.nseg;
  wseg = data.baseSeg.imgRegions;
  % compute foreground probability
  sumObj = zeros(nSps, 1); cntObj = zeros(nSps, 1);
  for k=1:numel(wseg)
    sumObj(wseg(k))=sumObj(wseg(k))+pm(k);
    cntObj(wseg(k))=cntObj(wseg(k))+1;
  end
  wObj = sumObj./cntObj; edgePairs = data.baseResult.chainLR;
  threshold = 0.4; 
  vWeight = (data.baseResult.edgeCost); 
  vWeight(min(wObj(edgePairs(:,1)), wObj(edgePairs(:,2)))<threshold) = 1e-8;
  vWeight(max(wObj(edgePairs(:,1)), wObj(edgePairs(:,2)))<threshold) = 1-1e-8;
  
  if opt.debug % merging weight map
    eim = zeros(size(data.depths));
    for k=1:numel(data.baseResult.edgeletChain)
      for l=data.baseResult.edgeletChain{k}
        eim(data.baseResult.edgeIndices{l}) = 1-vWeight(k);
      end
    end
    figure(2); clf; imshow(eim); drawnow;
  end
  
  wObj(wObj<threshold) = 0;
  t = wObj.*(cntObj/sum(cntObj)).^.5;
  t = t/sum(t);
  [~, seeds] = find(mnrnd(1, t, opt.maxNRegion*20));
  seeds = seeds(randperm(opt.maxNRegion*20));
  
  vWeight(vWeight<1e-8) = 1e-8;
  result.prop = mexMyProposal(int16(edgePairs), vWeight, cntObj/sum(cntObj), int16(seeds));
  isValid = true(size(result.prop));
  for ii=1:numel(result.prop)
    if sum(cntObj(result.prop{ii}+1))<opt.too_small, isValid(ii)=false; end
  end
  
  result.prop = {result.prop{isValid}}; seeds = seeds(isValid);
  isValid = true(size(result.prop));
  % take unique
  for ii=1:numel(result.prop)
    if ~isValid(ii), continue; end
    for jj=(ii+1):numel(result.prop)
      isValid(jj) = isValid(jj) & ~isequal(result.prop{ii}, result.prop{jj});
    end
  end
  result.prop = {result.prop{isValid}}; seeds = seeds(isValid);
  
  occ = false(nSps, numel(result.prop));
  areaObj = zeros(nSps,1);
  
  for jj=1:numel(result.prop)
    occ(result.prop{jj}+1, jj)=true;
    areaObj(jj)=sum(cntObj(result.prop{jj}+1));
  end
%  IOU = ones(numel(result.prop));
%   for ii=1:nSps
%     occSeg = find(occ(ii,:));
%     for j1=occSeg, for j2=occSeg, IOU(j1,j2)=IOU(j1,j2)+cntObj(ii); end; end
%   end
%   t = repmat(areaObj, [1 numel(result.prop)])+repmat(areaObj', [numel(result.prop) 1]);
%   IOU=IOU./(t-IOU);
  
  occ = bsxfun(@times, occ, cntObj);
  [IDX,C] = kmeans( occ', opt.maxNRegion, 'EmptyAction', 'singleton');
  [uIDX, ~, t]=unique2(IDX, 'legacy');
  C=C(uIDX, :);
  
  nProp = size(C, 1);
  ids = zeros(nProp, 1);
  for ii=1:nProp
    cl = find(IDX==uIDX(ii));
    [~, id] = min((sum(bsxfun(@minus, occ(:, cl), C(ii, :)').^2)));
    ids(ii) = cl(id);
  end
  result.prop = result.prop(ids);
  result.seeds = seeds(ids);
  result.seg = wseg;
  result.wObj = wObj;
  result.cntObj = cntObj;
  regionMasks = {};
  
  for j=1:nProp
    regionMasks{j} = ismember(result.seg, result.prop{j}+1);
    result.prop{j} = result.prop{j}+1;
  end
  data.regionProp = regionMasks;
  if opt.debug && 0
    figure(3); clf;
    for j=1:nProp
      I = data.images * .5 + .2; t = I(:,:,2); 
      m = regionMasks{j};
      t(wseg==result.seeds(j)) = t(wseg==result.seeds(j)) + .3; I(:,:,1)=t;
      t(m) = t(m) + .3; I(:,:,2)=t;
      imshow(I);
      pause;
    end
  end
  data.region_info = result;
  if opt.debug
    figure(3), clf; imagesc(sum(cat(3, regionMasks{:}), 3)); drawnow;
  end
end