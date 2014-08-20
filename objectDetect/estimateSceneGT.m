function data = estimateSceneGT(data, opt)

  % pre-compute min and max depth of each object
  prop = struct('minD', [], 'maxD', [], 'fvc', [], 'type', [], 'mask', [], 'dD', []); cnt = 0;
  result_layout = data.layoutProp.info;
  result_object = data.objectProp.info;
  itemGroup = []; unary = [];
  for ii=1:numel(result_layout)
    cnt = cnt + 1; fvc = result_layout(ii);
    minD = kinectMesa(fvc, false); maxD = kinectMesa(fvc, true);
    maxD = 10*(maxD>0); t = maxD-minD;
    prop(cnt) = struct('minD', minD, 'maxD', maxD, 'fvc', fvc, 'type', result_layout(ii).type, 'mask', ~~minD, 'dD', t);  
    itemGroup(cnt) = ii;
  end
  
  for ii=1:numel(result_object)
    if isempty(result_object(ii).retrieval), continue; end
    for jj=1:numel(result_object(ii).retrieval)
      fvc = struct('vertices', result_object(ii).retrieval(jj).vertices, 'faces', result_object(ii).retrieval(jj).faces);
      minD = kinectMesa(fvc); 
      maxD = kinectMesa(fvc, true);
      cnt = cnt + 1; t = maxD-minD;
      prop(cnt) = struct('minD', minD, 'maxD', maxD, 'fvc', fvc, 'type', result_object(ii).retrieval(jj).label, 'mask', ~~minD, 'dD', t);
      itemGroup(cnt) = ii+numel(result_layout);
    end        
  end

  regionMasks={result_object.mask};
  
  if opt.debug, fprintf('Computing Unary Potentials\n'); end
  % compute unary
  cnt = numel(result_layout);
  for ii=1:numel(result_object)
    if isempty(result_object(ii).retrieval), continue; end
    reg1mask = regionMasks{ii};
    for jj=1:numel(result_object(ii).retrieval)
      cnt = cnt+1;
      D = prop(cnt).minD;
      m = ~~D;
      unary(cnt) = ...
                sum(abs(data.depths(m&reg1mask)-D(m&reg1mask))) ...
              + sum(max(data.depths(m&~reg1mask)-D(m&~reg1mask),0)) ...
              + opt.missing_C * sum(sum(~m&reg1mask));
    end
  end
  
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
      overlap(ii,jj) = sum(m(:));
      overlap(jj,ii) = overlap(ii,jj);
    end
  end

  % setup linear program
  nNodes = numel(itemGroup);
  [n1, n2] = find(overlap>0);
  pairwise = overlap(overlap>0);
  nEdges = numel(n1);

  elem = unique(itemGroup);
  A1 = [];
  % nodes constraints
  cnt = 0;
  for ii=1:numel(elem)
    if ~ismember(ii, itemGroup), continue; end
    cnt = cnt + 1;
    A1(1:nNodes,cnt)=itemGroup==ii;
  end
  A1 = [A1; zeros(nEdges, cnt)];
  B1 = ones(cnt,1);
  % edge constraints
  A2 = zeros(nNodes+nEdges, nNodes);
  B2 = zeros(nNodes,1);
  for ii=1:nNodes
    A2(ii,ii)=-1;
    A2(nNodes+find(n1==ii),ii)=1;
  end
  % remove invalid constraints
  activeConst = (sum(A2)~=-1);
  A2 = A2(:, activeConst);
  B2 = B2(activeConst);
  LB = zeros(nNodes+nEdges,1);
  UB = ones(nNodes+nEdges,1);
  f = double(cat(2, unary, pairwise'));
  Aeq = cat(1, A1', A2'); Beq = cat(1, B1, B2);
  x = linprog(f,[],[],Aeq,Beq,LB,UB);
  x = x(1:nNodes);

  fvc = struct('vertices', [], 'faces', []);
  total_fvc = struct('vertices', [], 'faces', []);
  uid = unique(itemGroup);
  for ii=1:numel(uid)
    %[~, maxi] = max(x(itemGroup==ii));
    [~, maxi] = max(x'.*double(itemGroup==uid(ii)));
    fvc(ii).faces = prop(maxi).fvc.faces;
    fvc(ii).vertices = prop(maxi).fvc.vertices;
    fvc(ii).label = prop(maxi).type;

    total_fvc.faces = cat(1, total_fvc.faces, size(total_fvc.vertices,1)+prop(maxi).fvc.faces);
    total_fvc.vertices = cat(1, total_fvc.vertices, prop(maxi).fvc.vertices);
  end
  data.final_fvc = fvc;
  data.total_fvc = total_fvc;
end