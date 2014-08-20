% Extracts boundary features from the given imgRgbage.
%
% Region Features:
%   area - the number of pixels in the region (relative to the imgRgbage size).
% 
% Computes features for each edgelet based on boundary and geometry
% confidence imgRgbages.  This version (3) computes geometry confidences as
% means of superpixel values on either side of the edglet.
%
% Input:
%   boundaryInfo - structure of superpixel boundaries
%   imgRgb(H, W, 3) - color imgRgbage
%   pbImg(H, W, norient) - probability of boundary imgRgbage at each orientation
%   gconf(H, W, ngeoms) - confidence in each geometric label
%
% Output:
%   X(nboundaries, :) - features for the boundaries
%   Y(nboundaries, 2) - geometry on each size or 0 for no edge
function X = compile_boundary_classifier_feature_struct(boundaryInfo, info3d, imgRgb, pbImg)

  R = boundaryInfo.nseg;
  edges = boundaryInfo.edges;
  ne = boundaryInfo.ne;
  [H, W] = size(boundaryInfo.imgRegions);

  X.edge.pb = zeros(ne, 1);
  X.edge.theta = zeros(ne, 1);
  X.edge.thetaStart = zeros(ne*2, 1);
  X.edge.thetaEnd = zeros(ne*2, 1);
  X.edge.smoothness = zeros(ne, 1);
  X.edge.length = zeros(ne, 1);
  X.edge.chains = zeros(ne, 1); % chains;
  X.edge.edge2chain = zeros(ne, 1);
  X.edge.chainsize = zeros(ne, 1);

  X.region.colorMean = zeros(R, 3); % Lab color mean
  X.region.colorHist = [];  % set number of bins later (max 512)
  X.region.X = zeros(R, 3); % min, mean, max 3D position
  X.region.Y = zeros(R, 3);
  X.region.Z = zeros(R, 3);
  X.region.sample3D = zeros(R, 5*4); % random sample of 5 3D points (with appended one vector)
  X.region.norms3d = zeros(R, 10); % mean of abs normal values and 7 normals categories
  X.region.planes = zeros(R, max(info3d.planeMap(:)));
  X.region.planeparam = zeros(R, 4);
  X.region.x = zeros(R, 3); % min, mean, max 2D position
  X.region.y = zeros(R, 3);
  X.region.area = zeros(R, 1);

  %% Edge statistics

  % get discretize orientation into 1=1|2, 2=1/2, 3=1_2, 4=2\1
  theta = boundaryInfo.edges.thetaUndirected;
  rels = (theta < -3*pi/8) + (theta < -pi/8) + (theta < pi/8) +  (theta < 3*pi/8); 
  rels = mod(rels, 4) + 1;

  X.edge.theta = boundaryInfo.edges.thetaDirected;  % directed edge angle

  % pbImg(H, W, [ -, \, |, / ]) (direction of edge)
  pbmap = [3 4 1 2];
  for ii = 1 : size(pbImg, 3)
    pbImg(:,:,ii) = ordfilt2(pbImg(:,:,ii), 25, true(5));
  end

  % compute features
  for ii = 1 : ne
    eind = edges.indices{ii};

    X.edge.length(ii) = numel(eind); % edge length

    pbsubi = pbmap(rels(ii));    
    ind = double(eind + (pbsubi-1)*H*W);
    X.edge.pb(ii) = sum(pbImg(ind))/numel(ind);  % mean pb               

    % short-range angles
    y = mod(ind-1, H)+1;
    x = floor((ind-1)/H)+1;
    ni = numel(ind);
    de = 10; % length of edge used to measure angle
    x1 = x([1 min(de, ni)]);
    y1 = y([1 min(de, ni)]);
    x2 = x([max(ni-de+1, 1) ni]);
    y2 = y([max(ni-de+1, 1) ni]);
    X.edge.thetaStart(ii) = atan2(-(y1(2)-y1(1)), x1(2)-x1(1));
    X.edge.thetaEnd(ii) = atan2(-(y2(2)-y2(1)), x2(2)-x2(1));

    X.edge.smoothness(ii) = (numel(ind)-1) / (abs(x(end)-x(1)) + abs(y(end)-y(1))+1);    
  end

  X.edge.thetaStart(ne+1:end) = X.edge.thetaEnd(1:ne)+pi;
  X.edge.thetaEnd(ne+1:end) = X.edge.thetaStart(1:ne)+pi;

  thetaEnd = mod(X.edge.thetaEnd*180/pi, 360);
  thetaStart = mod(X.edge.thetaStart*180/pi, 360);

  % chain together edgelets 
  [chains, e2chain, chainsize] = chainEdgelets([X.edge.pb ; X.edge.pb], ...
      edges.adjacency, thetaStart, thetaEnd, 0.02, 45);
  X.edge.chains = chains;
  X.edge.edge2chain = single(e2chain);
  X.edge.chainsize = single(chainsize);

  %% Region statistics

  % get area and position stats
  stats = regionprops(boundaryInfo.imgRegions, 'PixelIdx', 'Area', 'Centroid', 'BoundingBox');

  % The number of pixels.
  area = cat(1, stats(:).Area);
  X.region.area = area / (H * W);

  bbox = cat(1, stats(:).BoundingBox);
  centroid = cat(1, stats(:).Centroid);
  minx = bbox(:,1);
  meanx = centroid(:, 1);
  maxx = minx + bbox(:,3);
  miny = bbox(:,2);
  meany = centroid(:, 2);
  maxy = miny + bbox(:,4);
  X.region.x = [minx meanx maxx]/W; % left, center, right
  X.region.y = (1-[maxy meany miny]/H); % bottom, center, top

  % get lab color imgRgbage
  imgRgb = RGB2Lab(imgRgb);

  % get discrete imgRgbage with nb values per channel
  nb = 8;
  mincols = repmat(min(min(imgRgb)), [H W]);
  maxcols = repmat(max(max(imgRgb))+1E-10, [H W]);
  imgRgbd = floor((imgRgb-mincols)./(maxcols-mincols)*nb);
  imgRgbd = imgRgbd(:, :, 1) + imgRgbd(:, :, 2)*nb + imgRgbd(:,:,3)*nb*nb + 1;

  % compute XYZ min/mean/max values and mean abs normal values and categories
  idx = {stats(:).PixelIdxList};
  XimgRgb = info3d.points3d(:, 1); 
  YimgRgb = info3d.points3d(:, 2);  
  ZimgRgb = info3d.points3d(:, 3); 
  absnorms3d = abs(info3d.normals);
  i1 = absnorms3d(:, 2)>0.95; % horizontal
  i2 = absnorms3d(:, 1).^2+absnorms3d(:, 3).^2>0.95.^2; % vertical
  i3 = absnorms3d(:, 1)>0.95; % aligned with wall 1
  i4 = absnorms3d(:, 3)>0.95; % aligned with wall 2
  i5 = i3 | i4; % aligned with either wall
  i6 = i2 &~i3 &~i4; % vertical but not aligned with any wall
  i7 = ~(i1 | i2 | i3); % not vertical or horizontal
  absnorms3d = [absnorms3d i1 i2 i3 i4 i5 i6 i7];
  for k = 1:R
      Xvals = XimgRgb(idx{k}); 
      X.region.X(k, :) = [min(Xvals) sum(Xvals)/area(k) max(Xvals)];
      Yvals = YimgRgb(idx{k}); 
      X.region.Y(k, :) = [min(Yvals) sum(Yvals)/area(k) max(Yvals)];
      Zvals = ZimgRgb(idx{k}); 
      X.region.Z(k, :) = [min(Zvals) sum(Zvals)/area(k) max(Zvals)];
      X.region.norms3d(k, :) = sum(absnorms3d(idx{k}, :))/area(k);    

      XYZ = [Xvals Yvals Zvals ones(size(Xvals))]; % get plane parameters: (aX + bY + cZ + d = 0)
      [eigv, l] = eig(XYZ'*XYZ);
      X.region.planeparam(k, :) = eigv(:, 1)' ./ sqrt(sum(eigv(1:3, 1).^2));

      rp = randperm(numel(Xvals)); 
      npts = size(X.region.sample3D, 2)/4;
      X.region.sample3D(k, :) = reshape(XYZ(rp(1:npts), :), [1 npts*4]); % random set of 3D points
  end

  % percent of pixels within each plane label
  for k = 1:R
    vals = info3d.planeMap(idx{k});
    for p = 1:size(X.region.planes, 2)
      X.region.planes(k, p) = sum(vals==p) / area(k);
    end
  end                

  % compute mean color
  imgRgb = reshape(imgRgb, [H*W 3]);
  for ii = 1 : 3
    cimgRgb = imgRgb(:,ii);    
    for jj = 1 : R
      X.region.colorMean(jj,ii) = sum(cimgRgb(idx{jj})) / area(jj);
    end
  end

  % compute histograms of color and texture
  imgRgbd =reshape(imgRgbd, [H*W 1]);
  imgRegions = boundaryInfo.imgRegions;
  colorHist = zeros(R, nb*nb*nb);
  for k = 1:H*W
    s = imgRegions(k);
    colorHist(s, imgRgbd(k)) = colorHist(s, imgRgbd(k)) + 1;
  end

  keep = sum(colorHist, 1)>0;  % only keep bins that have at least one value
  colorHist = colorHist(:, keep);
  X.region.colorHist = single(colorHist ./ repmat(area, [1 sum(keep)]));
end
