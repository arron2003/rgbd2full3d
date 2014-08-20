function data = estimateScenePROP(data, opt) 
    % pre-compute min and max depth of each object
  prop = struct('minD', [], 'maxD', [], 'fvc', [], 'type', [], 'reg', [], 'vol', [], 'mask', [], 'dD', []); 
  cnt = 0;
  result_layout = data.layoutProp.info;
  result_object = data.objectProp.info;
  unary = []; nLayout = numel(result_layout);
  for ii=1:numel(result_layout)
    cnt = cnt + 1; fvc = result_layout(ii);
    minD = single(kinectMesa(fvc, false)); maxD = single(kinectMesa(fvc, true));
    maxD = 10*(maxD>0); t = maxD-minD;
    prop(cnt) = struct('minD', minD, 'maxD', maxD, 'fvc', fvc, 'type', result_layout(ii).type,...
      'reg', 0, 'vol', sum(t(:)), 'mask', ~~minD, 'dD', t);  
  end
  
  for ii=1:numel(result_object)
    if isempty(result_object(ii).retrieval), continue; end
    for jj=1:numel(result_object(ii).retrieval)
      fvc = struct('vertices', result_object(ii).retrieval(jj).vertices, 'faces', result_object(ii).retrieval(jj).faces);
      minD = kinectMesa(fvc); 
      maxD = kinectMesa(fvc, true);
      cnt = cnt + 1; t = maxD-minD;
      prop(cnt) = struct('minD', minD, 'maxD', maxD, 'fvc', fvc, 'type', result_object(ii).retrieval(jj).label,...
        'reg', ii, 'vol', sum(t(:)), 'mask', ~~minD, 'dD', t);
    end        
  end

  regionMasks={result_object.mask};
  
  if opt.debug, fprintf('Computing Unary Potentials\n'); end
  % compute unary
  cnt = 0; unary = zeros(1, numel(prop));
  layoutP = 1-data.pm(:,:,4);
  objP = data.pm(:,:,4);
  md = max(data.depths(:));
  for ii=1:numel(result_layout)
    cnt = cnt+1;
    D = prop(cnt).minD; m = ~~D;
    t = (m & abs(D-data.depths)<opt.nms).*layoutP;
    unary(cnt) = - sum(t(:)) - opt.nms*sum(md-D(m));
  end
  
  cnt = numel(result_layout);
  for ii=1:numel(result_object)
    if isempty(result_object(ii).retrieval), continue; end
    reg1mask = regionMasks{ii};
    for jj=1:numel(result_object(ii).retrieval)
      cnt = cnt+1;
      D = prop(cnt).minD;
      m = ~~D;
      t = (abs(data.depths-prop(cnt).minD)<opt.nms).*objP;
      unary(cnt) = -opt.nms*sum(t(:)) - opt.nms*prop(cnt).vol...
              + sum(max(data.depths(m)-D(m), 0));
    end
  end

  itemGroup = [prop.reg];
  itemGroup(1:nLayout) = 1:nLayout;
  itemGroup((nLayout+1):end) = itemGroup((nLayout+1):end)+nLayout;
  if opt.debug, fprintf('Computing Pairwise Potentials\n'); end
  
  % compute pairwise
  
  overlap = zeros(numel(itemGroup));
  for ii=1:numel(itemGroup)
    for jj=(ii+1):numel(itemGroup)
      if jj<=numel(result_layout), continue; end % only consider object-* overlap
      t = prop(ii).mask & prop(jj).mask;
      if (~any(t(:))), continue; end
      m = zeros(size(prop(jj).maxD));
      m(t) = min( prop(jj).maxD(t) - prop(ii).minD(t), prop(ii).maxD(t) - prop(jj).minD(t) );
      m(t) = min(m(t), prop(ii).dD(t));
      m(t) = min(m(t), prop(jj).dD(t));
      m(t) = max(m(t), 0);
      if ii<=nLayout 
        m = max(0, m-opt.nms);
        overlap(ii,jj) = sum(m(:));
      else
        overlap(ii,jj) = sum(m(:));
        racing=abs(prop(ii).minD(t)-prop(jj).minD(t));
        overlap(ii,jj) = overlap(ii,jj) + sum(racing<opt.nms);
      end
      
      % penalize racing condition
      if prop(jj).reg==prop(ii).reg % no objects can be from the same region proposals
        overlap(ii,jj) = 1e5;
      end
      
      overlap(jj,ii) = overlap(ii,jj);
    end
  end

  nProp = numel(prop);
  nObject = numel(prop) - nLayout;
  ios = overlap;
  for i=1:nProp
    for j=1:nProp
      ios(i,j) = ios(i,j) / min(prop(i).vol, prop(j).vol);
    end
  end
  ios(ios>1e3)=0;
  ios(isnan(ios))=0;

  % completeness constraints
  wseg = data.baseSeg.imgRegions;
  nSeg = data.baseSeg.nseg;
  cntObj = data.region_info.cntObj;
  wObj = data.region_info.wObj;
  validRegion = find(cntObj>opt.too_small & wObj>0.4);
  C1 = zeros(numel(prop), numel(cntObj));
  for jj=(numel(result_layout)+1):numel(prop)
    m = abs(prop(jj).minD-data.depths)<opt.nms;
    t0 = histc(wseg(m), 1:nSeg);
    t0 = t0(:).*wObj(:);
    t = find(t0); t0=t0(t);
    [t1, t] = ismember(t, validRegion);
    t1 = t0(t1);
    t = t(~~t);
    C1(jj, t) = t1;
  end

  %keyboard;
  [x, f] = beamSearch(unary, overlap, zeros(size(C1)), ios>.2);
  final_fvc = struct('vertices', [], 'faces', []); cnt = 0;
  for ii=find(x>.5)'
    cnt = cnt+1;
    final_fvc(cnt) = struct('vertices', prop(ii).fvc.vertices, 'faces', prop(ii).fvc.faces);
  end
  data.final_fvc = final_fvc;
  %keyboard;
end
