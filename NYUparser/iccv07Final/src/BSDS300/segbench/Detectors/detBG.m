function [bg,theta] = detBG(im,radius,norient)
% function [bg,theta] = detBG(im,radius,norient)
%
% Compute smoothed but not thinned BG fields.

if nargin<2, radius=0.01; end
if nargin<3, norient=8; end

[h,w,unused] = size(im);
idiag = norm([h w]);
if size(im, 3)==3, im=rgb2gray(im); end

% compute brightness gradient
[bg,theta] = cgmo(im,idiag*radius,norient,...
                  'smooth','savgol','sigmaSmo',idiag*radius);


