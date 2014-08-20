% Creates a visualization of structure classes.
%
% Args:
%   imgRgb - HxWx3 RGB image.
%   imgRegions - HxW region map, where each pixel ranges from 0 (missing
%                region) to R, the number of total regions.
%   structureLabelsGt - Rx1 vector of ground truth structure class labels.
%   structureLabelsPred - Rx1 vector of predicted structure class labels.
function vis_structure_classes(imgRgb, imgRegions, ...
    structureLabelsGt, structureLabelsPred)

  imgRegions = double(imgRegions);
  
  % colors for floor, support, objects, structures.
  hues = [0.83, 0.17, 0.75, 0.56];
  
  % adjusts brightness of colors.
  satlevel = 0.75;

  [H, W] = size(imgRegions);
  R = max(imgRegions(:));
  [xim, yim] = meshgrid(1:W, 1:H);

  % assign colors based on metalabels
  hueim = zeros(size(imgRegions));
  satim = zeros(size(imgRegions));
  labim = zeros(size(imgRegions));
  
  % Create the stripe mask.
  imgStripes = ones(H, W);
  stride = 10;
  stripeTemplate = ones(1,W);
  
  thickness = 5;
  for tt = 1 : thickness
    stripeTemplate(tt:stride:W) = 0;
  end
  
  for hh = 1 : H
    shift = mod(hh-1, stride)+1;
   imgStripes(hh,:) = circshift(stripeTemplate', shift)';
  end
%   imgStripes = ~imgStripes;
  imgStripes = logical(imgStripes);
  
  for rr = 1 : R
    regionMask = imgRegions == rr;
    if structureLabelsPred(rr) ~= structureLabelsGt(rr)
      regionMask = logical(regionMask .* imgStripes);
    end
    labim(regionMask) = structureLabelsPred(rr);
  end
  
  candycane = true; thickness = 1; spacing = 5; % whether to use candycane pattern, thickness and spacing
  
%   labim(imgRegions>0) = mlabels(imgRegions(imgRegions>0));
  for y = 1:4
    hueim(labim==y) = hues(y);
    satim(labim==y) = satlevel;
  end
  dispim = hsv2rgb(cat(3, hueim, satim, im2double(rgb2gray(imgRgb))*0.5+0.5));

  % draw boundaries
  [gx, gy] = gradient(imgRegions);
  e = abs(gx)+abs(gy)>0;
  dispim(repmat(e, [1 1 3])) = 1;
  if candycane  
    tmp = e & (mod(xim+yim, spacing)<thickness); 
    dispim(:, :, 2:3) = dispim(:, :, 2:3).*repmat(tmp==0, [1 1 2]);
  end
  
  imagesc(dispim);
end