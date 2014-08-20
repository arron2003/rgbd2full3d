% Gets the vector of labels assigned to each region in imgRegions.
%
% Args:
%   imgRegions - HxW image whose values range from 0 (ignored region) to R
%                (the maximum region value).
%   imgObjectLabels - HxW image whose values range from 0 (missing label)
%                     to C (the maximum label).
%   imgInstanceLabels - (optional) HxW image whose values range from 0 ...
%                       (missing label) to N (the maximum instance of that
%                       class).
%
% Returns:
%   objectLabels - a list of object class labels (1..C) for each region
%                  in the image.
%   instanceLabels - a list of each instance label (1..N) for each region
%                    in the image.
%   intersection - the percentage that each region overlaps with each
%                  object instance.
function [objectLabels, instanceLabels, intersection] = ...
    get_labels_from_instances(imgRegions, imgObjectLabels, imgInstanceLabels)
  
  assert(all(size(imgRegions) == size(imgObjectLabels)));
  assert(all(size(imgRegions) == size(imgInstanceLabels)));

  R = max(imgRegions(:));

  [instanceMasks, instanceObjectLabels] = get_instance_masks(...
    imgObjectLabels, imgInstanceLabels);
  N = numel(instanceObjectLabels);

  imgInstances = zeros(size(imgRegions));
  for ii = 1 : N
    imgInstances(instanceMasks(:,:,ii)) = ii;
  end

  imgInstances(imgObjectLabels == 0) = N + 1;

  [~, intersection, ~] = get_region_overlap(imgRegions, imgInstances);

  [~, inds] = max(intersection, [], 2);
  hasLabel = inds < N + 1;

  objectLabels = zeros(R, 1);
  objectLabels(hasLabel) = instanceObjectLabels(inds(hasLabel));
  instanceLabels = inds;

  numPixels = get_num_pixels_per_region(imgRegions);
  intersection = intersection ./ repmat(numPixels, [1, N+1]);
end
