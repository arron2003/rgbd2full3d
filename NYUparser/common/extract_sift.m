% Extracts sift descriptors from the given coordinates.
%
% This code is a very slight adaptation of the SIFT code from Svetlana
% Lazebnik's Spatial Pyramid code:
%   http://www.cs.illinois.edu/homes/slazebni/
%   http://www.cs.illinois.edu/homes/slazebni/research/SpatialPyramid.zip
%
% Args:
%   img - the input image, either a grayscale or rgb image. If the image is
%       RGB, its average over the 3 channels will be used. Should be a
%       double between 0 and 1.
%   coords - Nx2 set of (X,Y) coordinates.
%   params - parameters for SIFT.
%
% Returns:
%   sift_arr - an NxD (Nx128) matrix.
%   Z - the normalizing coefficient.
function [sift_arr, Z] = extract_sift(img, coords, params)
  error(nargchk(3,3,nargin));
  assert(isa(img, 'double') || isa(img, 'single'));
  N = size(coords, 1);
  
  %%%%%%%%%%%%%%
  % PATCH SIZE %
  %%%%%%%%%%%%%%
  patchSize = params.patchSize;
  assert(patchSize > 0, 'Patch Size must be positive');
  
  %%%%%%%%%%%%%%
  % NUM ANGLES %
  %%%%%%%%%%%%%%
  if isfield(params, 'numOris')
    assert(params.numOris > 1, 'Must have more than one orientation');
    numOris = params.numOris;
  else
    numOris = 8;
  end
  
  %%%%%%%%%%%%
  % NUM BINS %
  %%%%%%%%%%%%
  if isfield(params, 'numBins')
    assert(params.numBins > 0, 'Number of bins must be positive')
    numBins = params.numBins;
  else
    numBins = 4;
  end
  totalBins = numBins * numBins;
  
  %%%%%%%%%%%%%%%%%
  % Normalization %
  %%%%%%%%%%%%%%%%%
  if isfield(params, 'normMethod')
    normMethod = params.normMethod;
  else
    normMethod = 1;
  end
  
  alpha = 9; %% parameter for attenuation of angles (must be odd)
  
  img = mean(img,3);
  img = img / max(img(:));

  if nargin < 5
    sigma_edge = 1;
  end

  angle_step = 2 * pi / numOris;
  angles = 0:angle_step:2*pi;
  angles(numOris+1) = []; % bin centers

  [hgt wid] = size(img);

  [GX, GY] = gen_dgauss(sigma_edge);

  % add boundary (mirror):
  img = [img(2:-1:1,:,:); img; img(end:-1:end-1,:,:)];
  img = [img(:,2:-1:1,:) img img(:,end:-1:end-1,:)];

  img = img - mean(img(:));
  img_X = filter2(GX, img, 'same'); % vertical edges
  img_Y = filter2(GY, img, 'same'); % horizontal edges

  img_X = img_X(3:end-2,3:end-2,:);
  img_Y = img_Y(3:end-2,3:end-2,:);

  img_mag = sqrt(img_X.^2 + img_Y.^2); % gradient magnitude
  imgTheta = atan2(img_Y,img_X);
  imgTheta(isnan(imgTheta)) = 0; % necessary????

  % make orientation images
  imgOris = zeros([hgt, wid, numOris], 'single');

  % for each histogram angle
  imgCos = cos(imgTheta);
  imgSin = sin(imgTheta);
  for a=1:numOris
      % compute each orientation channel
      tmp = (imgCos*cos(angles(a))+imgSin*sin(angles(a))).^alpha;
      tmp = tmp .* (tmp > 0);

      % weight by magnitude
      imgOris(:,:,a) = tmp .* img_mag;
  end

  % Convolution formulation:
  r = patchSize/2;
  cx = r - 0.5;
  sample_res = patchSize/numBins;
  weight_x = abs((1:patchSize) - cx)/sample_res;
  weight_x = (1 - weight_x) .* (weight_x <= 1);

  for ii = 1 : numOris
    imgOris(:,:,ii) = conv2(weight_x, weight_x', imgOris(:,:,ii), 'same');
  end

  quarterBin = patchSize / 4;
  [sample_x, sample_y] = meshgrid(-quarterBin*1.5 : quarterBin : quarterBin * 1.5);
  sample_x = sample_x(:);
  sample_y = sample_y(:);
  
  sift_arr = zeros([N numOris * totalBins], 'single');
  b = 0;
  for jj = 1 : totalBins
    for nn = 1 : N
      sift_arr(nn, b+1:b+numOris) = imgOris(round(coords(nn,2)+sample_y(jj)), round(coords(nn,1)+sample_x(jj)), :);
    end
    b = b + numOris;
  end
  clear imgOris
  
  switch normMethod
    case 1
      % normalize SIFT descriptors
      ct = .1;
      sift_arr = sift_arr + ct;
      Z = sqrt(sum(sift_arr.^2, 2));
      sift_arr = sift_arr ./ repmat(Z, [1 size(sift_arr,2)]);

      Z = Z(:);
    case 2
      [sift_arr, Z] = sp_normalize_sift(sift_arr);
  end
end

function [GX, GY] = gen_dgauss(sigma)
  % laplacian of size sigma
  G = gen_gauss(sigma);
  [GX,GY] = gradient(G); 

  GX = GX * 2 ./ sum(sum(abs(GX)));
  GY = GY * 2 ./ sum(sum(abs(GY)));
end

function G = gen_gauss(sigma)
  if all(size(sigma)==[1, 1])
    % isotropic gaussian
    f_wid = 4 * ceil(sigma) + 1;
    G = fspecial('gaussian', f_wid, sigma);
  else
    % anisotropic gaussian
    f_wid_x = 2 * ceil(sigma(1)) + 1;
    f_wid_y = 2 * ceil(sigma(2)) + 1;
    G_x = normpdf(-f_wid_x:f_wid_x,0,sigma(1));
    G_y = normpdf(-f_wid_y:f_wid_y,0,sigma(2));
    G = G_y' * G_x;
  end
end
