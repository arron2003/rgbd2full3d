% Infers the assignment of points to planes using graph cuts. Note that the
% 3D point cloud in this function considers the XZ plane to be the ground
% plane.
%
% Args: 
%   imgRgb - HxWx3 double-valued matrix representing the RGB image.
%   planeNdxs - cell array of the indices assigned to each plane.
%   points3d - an Nx4 matrix, the 3d point cloud in homogeneous
%              coordinates.
%   normals - an Nx3 matrix, the surface normals at each point of the
%             image.
%   imgDepth - HxW depth image.
%   depthWeight - the amount to weight the depth unary term.
%
% Returns:
%   labelMap - HxW map of plane labels running from 0 (no plane) to P where
%              P is the total number of planes.
%   planes - 
function [imgLabels, planes] = cut_planes_registered(imgRgb, planeNdxs, ...
    points3d, normals, imgDepth, depthWeight)

  DEBUG = false;
  
  % parameters
  gaussBlur = fspecial('gaussian', 9, 1);
  settings.lambda1 = 1; %0.5;  %Ising Prior weight (5 by default)
  settings.lambda2 = 45; %45; %Edge sensitive weight (90 by defaul)
  unary_bound = 50;
  colorWeight = 0;
  distthresh = 0.005; % above this distance, do not allow depth to favor plane

  settings.GMM_nComponents = 5; %5; %Number of GMM components
  settings.GrabCut_iterations = 5; %5; %Number of times the color models are updated
  settings.neighboorSystem = 1;
  settings.p_norm = 1; % norm - note: it can also have multiple values. Example: settings.p = [1 2 3];
  settings.width = 1; %Should be an odd number. it imposes a square width*width
  options.nComponents = settings.GMM_nComponents; 
  options.iterations = 4;
  options.init = 0; %Initialisation Option: 0 - smart, 1 - random
  options.method = 'kmeans';                      				    
  options.covarType = 'full'; %'spherical';%
  options.min_var = (1/255)^2.0;
  options.LazySnapping = 0; 

  % initialize
  [H, W, ~] = size(imgRgb);
  P = numel(planeNdxs);

  imgRgb = imfilter(imgRgb, gaussBlur);
  imgBeta = Compute_beta(imgRgb,settings.neighboorSystem);

  % get initial label map in imgRgb space
  hasbeenset = false(H, W);
  imgLabels = zeros(H, W);
  for pp = 1 : P
    tmp = false(H, W);
    tmp(planeNdxs{pp}) = true;
    
    % Set pixels that have multiple labels to background
    imgLabels(tmp & hasbeenset) = 0;
    imgLabels(tmp & ~hasbeenset) = pp;
    %planeNdxs{p} = planeNdxs{p}(~hasbeenset(planeNdxs{p}));
    hasbeenset(tmp) = true;    
  end
  
  imgLabels(imgLabels == 0) = P+1;
  
  if DEBUG
    figure;
    imagesc(imgLabels);
    axis image;
    colormap jet;
    title('Initial Label Map [hit any key to continue]');
    pause;
  end

  % estimate color distributions
  [logpcolor_bg, logpcolor_fg] = estimate_color_dist(H, W, P, ...
    imgRgb, imgLabels, colorWeight, options);

  % estimate distributions of distance to plane vs bg.
  [logpdist_bg, logpdist_fg, planes, planesN] = ...
    estimate_dist_to_plane_distrib(H, W, P, planeNdxs, points3d, imgDepth, distthresh);


  % estimate surface normal distributions   
  [logpsurf_bg, logpsurf_fg] = estimate_surf_normal_dists(...
    H, W, P, planeNdxs, normals, planesN);

  % compute unaries, mapping unaries on 3d points to the image
  unary_fg = colorWeight*clipPotentials(-logpcolor_fg + logpcolor_bg, unary_bound);        
  for pp = 1 : P
    distpot = (-logpdist_fg(:,pp)+logpdist_bg(:,pp)).*depthWeight(:);
    unary_fg(:,:,pp) = unary_fg(:,:,pp) + clipPotentials(reshape(distpot, [H W]), unary_bound*10*depthWeight);
    surfpot = (-logpsurf_fg(:,pp)+logpsurf_bg(:,pp)).*depthWeight(:);
    unary_fg(:,:,pp) = unary_fg(:,:,pp) + clipPotentials(reshape(surfpot, [H W]), unary_bound*2*depthWeight);        
  end

  iter = 0;
  changed = true;
  while changed

    iter = iter + 1;
    lastmap = imgLabels;

    % perform a set of alpha expansions
    ordering = randperm(P+1);
    for alpha = ordering
      unary_alpha = zeros(H, W, 2);        
      for pp = 1 : P + 1
        if pp == alpha
          unary_alpha(imgLabels == pp) = -100000;
        else                
          ind = find(imgLabels == pp);
          indp = H * W * (pp-1) + ind;
          inda = H * W * (alpha-1) + ind;
          unary_alpha(imgLabels == pp) = unary_fg(inda)-unary_fg(indp);
        end            
      end    

      gcseg = mex_graphCut(unary_alpha, settings.lambda1, ...
          settings.lambda2, settings.neighboorSystem, imgRgb, imgBeta);     
      imgLabels(gcseg == 0) = alpha;
    end

%     if 0 
%         tmpunary = reshape(unary_fg, [H*W P+1]);
%         tmp = 0; for k = 1:numel(imgLabels), tmp = tmp + tmpunary(k, imgLabels(k)); end
%         fprintf('Avg Unary: %0.3f\n', tmp/H/W);
%     end

    if DEBUG
      figure(100);
      clf;
      
      subplot(1,2,1);
      imagesc(imgRgb);
      axis image
      
      subplot(1,2,2);
      imagesc(imgLabels);
      axis image;
      colormap jet;
      title(sprintf('Label Map iteration %d', iter));
      pause(.3);
    end

    % get new label indices of 3d points.
    for pp = 1 : P
      planeNdxs{pp} = find(imgLabels == pp);
    end

    changed = any(lastmap(:) ~= imgLabels(:));
  end

  if DEBUG
    figure;
    imagesc(imgLabels);
    title('Final image map [press any key to continue]');
    pause;
  end
end

function [logpcolor_bg, logpcolor_fg] = estimate_color_dist(H, W, P, ...
  imgRgb, imgLabels, colorWeight, options)

  % estimate color distributions
  logpcolor_bg = zeros(H, W, P+1);
  logpcolor_fg = zeros(H, W, P+1);
  if colorWeight > 0
    for pp = 1 : P+1
      imgRgbvalues = reshape(imgRgb, [H*W 3]);
      [pF, mF, CF] = ComputeGMM_cr(imgRgbvalues(imgLabels == pp, :)', options, 1, 0, 0, 0, 0);
      [pB, mB, CB] = ComputeGMM_cr(imgRgbvalues(imgLabels ~= pp, :)', options, 1, 0, 0, 0, 0);                
      logpcolor_fg(:,:,pp) = -reshape(EvGMM_sara(imgRgbvalues', pF, mF, CF), [H W])/2;
      logpcolor_bg(:,:,pp) = -reshape(EvGMM_sara(imgRgbvalues', pB, mB, CB), [H W])/2;
    end    
  end
end

function [logpdist_bg, logpdist_fg, planes, planesN] = ...
    estimate_dist_to_plane_distrib(H, W, P, planeNdxs, points3d, ...
      imgDepth, distthresh)

  % estimate distributions of distance to plane vs bg   
  planes = zeros(P, 4);
  planesN = zeros(P, 3);         
  logpdist_bg = zeros(H*W, P);
  logpdist_fg = zeros(H*W, P);
  for p = 1 : P

    A = points3d(planeNdxs{p}, :);
    
    % Do NOT remove the empty argument here. Matlab returns just the
    % eigenvalues instead of the eigenvectors otherwise.
    [eigv, ~] = eig(A'*A);        
    planes(p, :) = eigv(:, 1)';
    planes(p, :) = planes(p, :) / sqrt(sum(planes(p, 1:3).^2));
    planesN(p, :) = planes(p, 1:3); %./sqrt(sum(planes(p, 1:3).^2));                        

    dist = points3d*planes(p, :)';
    dist = abs(dist) ./ (imgDepth(:)+0.05); %.^1.5;

    dx =  0.0005:0.001:0.0995;
    hf = hist(dist(planeNdxs{p}), dx); 
    for k = 2:numel(hf), hf(k) = min(hf(k), hf(k-1)); end;
    hf = hf + 1; hf = hf / sum(hf);                         

    bin = min(max(ceil(dist*1000),1), 100);

    bgpts = true(H*W, 1);
    bgpts(planeNdxs{p}) = false;         
    hb = hist(dist(bgpts), dx); 
    hb = hb + 1; hb = hb/sum(hb);                                
    ind = dx > distthresh;
    hf(ind) = min(hf(ind), hb(ind)); % disallow large distances        

    logpdist_fg(:, p) = log(hf(bin));
    logpdist_bg(:, p) = log(hb(bin));
    %figure(1), hold off, plot(dx, log(hf)-log(hb)); pause;    
  end
end


function [logpsurf_bg, logpsurf_fg] = estimate_surf_normal_dists(...
    H, W, numPlanes, planeNdxs, normals, planesN)
  
  logpsurf_bg = zeros(H*W, numPlanes);
  logpsurf_fg = zeros(H*W, numPlanes);
  for ii = 1 : numPlanes
    Ndist = 1 - abs(normals*planesN(ii, :)');
    dx =  0.005 : 0.01 : 0.995 ;
    bin = max(ceil(Ndist*100),1);               
    hf = hist(Ndist(planeNdxs{ii}), dx); 
    for k = 2 : numel(hf)
      hf(k) = min(hf(k), hf(k-1));
    end
    hf = hf + 1;        
    hf = hf / sum(hf); 
    logpsurf_fg(:, ii) = log(hf(bin));

    bgpts = true(H*W, 1);
    bgpts(planeNdxs{ii}) = false;
    hb = hist(Ndist(bgpts), dx); 
    hb = hb + 1;
    hb = hb / sum(hb); 
    logpsurf_bg(:, ii) = log(hb(bin));
    
%     if DEBUG
%       figure(1);
%       hold off;
%       plot(dx, log(hf)-log(hb));
%       title(sprintf('Plane %d/%d', ii, numPlanes));
%       pause;
%     end
  end
end

function pot = clipPotentials(pot, bound)
  pot = max(min(pot, bound), -bound);
end
