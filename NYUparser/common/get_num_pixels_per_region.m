% Counts the number of pixels in each region.
%
% Args:
%   imgRegions - HxW image whose values range from 0 (ignored region) to R
%                (the maximum region value).
%
% Returns:
%   numPixels - Rx1 vector listing the number of pixels in each region.
function numPixels = get_num_pixels_per_region(imgRegions)
  R = max(imgRegions(:));
  numPixels = mex_num_pixels_per_region(int32(imgRegions), R);
end
