% Args:
%   scene - the name of the image's scene.
%   imageNdx - the frame number
%   imgRegions - a map of regions. Any pixels labeled as 0 will be ignored.
function features = get_surf_normal_features(regionMasks, points3d, normals, params)
  Consts;

  R = size(regionMasks,3);
  
  % Get the normals for the scene in spherical coordinates. Thetas are
  % between 0 and pi
  [inclinations, angles] = get_spherical_coord_values(points3d, normals);
  
  % Setup the output.
  numSpatialBins = get_spatial_bin_count_1d(params.numLevels);
  numSpatialBins = numSpatialBins * 2;
  D = (params.numInclinationBins + params.numAngleBins) * numSpatialBins;
  features = zeros(R, D);
  
  % Setup the bins.
  inclinationBins = 0  : pi/params.numInclinationBins : pi;
  inclinationBins(end) = inclinationBins(end) + 1; % Hack for histc.
  
  angleBins = 0 : 2*pi/params.numAngleBins : 2*pi;
  angleBins(end) = angleBins(end) + 1; % Hack for histc
  
  for rr = 1 : R
    regionMask = regionMasks(:,:,rr);
    
    % Bin the surface normals in terms of heights first.
    regionPoints3d = points3d(regionMask, :);
    featuresHeight = get_height_pyramid(regionPoints3d, inclinations(regionMask), inclinationBins, ...
        angles(regionMask), angleBins, params.numLevels);
    assert(~any(isnan(featuresHeight)));
      
    % Now, bin the surface normals across the principle axis in the XY
    % plane.
    regionPoints3d = points3d(regionMask, :);
    featuresLength = get_length_pyramid(regionPoints3d, inclinations(regionMask), inclinationBins, ...
        angles(regionMask), angleBins, params.numLevels);  
    assert(~any(isnan(featuresLength)));
    
    features(rr,:) = [featuresHeight featuresLength];
  end
end

% Args:
%   regionPoints3d - Mx3 point cloud.
%   inclinations - Mx1 inclinations.
%   incBins - inclination bins.
%   angles - Mx1 angles in the XY plane.
%   angBins - angle bins.
%   numLevels - the number of levels in the pyramid.
function features = get_height_pyramid(regionPoints3d, inclinations, incBins, ...
    angles, angBins, numLevels)
  
  % Initialize the features
  numSpatialBins = get_spatial_bin_count_1d(numLevels);
  H = numel(incBins) + numel(angBins) - 2;
  D = H * numSpatialBins; 
  features = zeros(1, D);
  
  maxHeight = max(regionPoints3d(:,3));
  minHeight = min(regionPoints3d(:,3));
  regionHeight = maxHeight - minHeight;
  
  % Pointers to the current location in the feature array.
  fp = 1;
  
  for level = 1 : numLevels
    numBins = 2^(level-1);
    heightOfBin = regionHeight / numBins;
    
    for bin = 1 : numBins
      startHeight = minHeight + (bin-1) * heightOfBin;
      endHeight = minHeight + bin * heightOfBin;
      
      % TODO: sensitive to outliers
      inds = regionPoints3d(:,3) >= startHeight & regionPoints3d(:,3) < endHeight;
      if nnz(inds) > 0
      
        incHist = histc(inclinations(inds), incBins);
        assert(incHist(end) == 0);
        incHist(end) = [];
        incHist = incHist ./ sum(incHist);

        angHist = histc(angles(inds), angBins);
        assert(angHist(end) == 0);
        angHist(end) = [];
        angHist = angHist ./ sum(angHist);

        features(fp : fp + H - 1) = [incHist(:)' angHist(:)'];
      end
      
      fp = fp + H;
    end
  end
end

% Args:
%   regionPoints3d - Mx3 point cloud.
%   inclinations - Mx1 inclinations.
%   incBins - inclination bins.
%   angles - Mx1 angles in the XY plane.
%   angBins - angle bins.
%   numLevels - the number of levels in the pyramid.
function features = get_length_pyramid(regionPoints3d, inclinations, incBins, ...
    angles, angBins, numLevels)
  
  % Initialize the features
  numSpatialBins = get_spatial_bin_count_1d(numLevels);
  H = numel(incBins) + numel(angBins) - 2;
  D = H * numSpatialBins; 
  features = zeros(1, D);
  
  % Find the principle axis in XY
  prinAxis = princomp(regionPoints3d(:,1:2));
  prinAxis = prinAxis(:,1);
  
  % Project all points onto this axis.
  N = size(regionPoints3d, 1);
  coeffs = (regionPoints3d(:,1:2) - repmat(mean(regionPoints3d(:,1:2)), [N 1])) * prinAxis;
  
  maxCoeff = max(coeffs);
  minCoeff = min(coeffs);
  regionLength = maxCoeff - minCoeff;
  
  % Pointers to the current location in the feature array.
  fp = 1;
  
  for level = 1 : numLevels
    numBins = 2^(level-1);
    lengthOfBin = regionLength / numBins;
    
    for bin = 1 : numBins
      startCoeff = minCoeff + (bin-1) * lengthOfBin;
      endCoeff = minCoeff + bin * lengthOfBin;
      
      inds = coeffs >= startCoeff & startCoeff < endCoeff;
      if nnz(inds) > 0
        incHist = histc(inclinations(inds), incBins);
        assert(incHist(end) == 0);
        incHist(end) = [];
        incHist = incHist ./ sum(incHist);

        angHist = histc(angles(inds), angBins);
        assert(angHist(end) == 0);
        angHist(end) = [];
        angHist = angHist ./ sum(angHist);

        features(fp : fp + H - 1) = [incHist(:)' angHist(:)'];
      end
      fp = fp + H;
    end
  end
end

function numSpatialBins = get_spatial_bin_count_1d(numLevels)
  assert(numLevels > 0, 'numLevels must be positive');
  numSpatialBins = 0;
  for ii = 1 : numLevels
    numSpatialBins = numSpatialBins + 2^(ii-1);
  end
end

function normals = flip_normals_towards_viewer(normals, points)
  points = points ./ repmat(sqrt(sum(points.^2, 2)), [1, 3]);
  
  proj = sum(points .* normals, 2);
  
  flip = proj > 0;
  normals(flip, :) = -normals(flip, :);
end


function [imgThetas, imgPsis] = get_spherical_coord_values(points3d, normals) 
  Consts;

  [~, sz] = get_projection_mask();
    
  normals = flip_normals_towards_viewer(normals, points3d);
  
  [~, thetas, psis] = cart2sphere(normals(:,1), normals(:,2), normals(:,3));
  thetas(isnan(thetas)) = 0;
  
  imgThetas = reshape(thetas, sz);
  imgPsis = reshape(psis, sz);
end
