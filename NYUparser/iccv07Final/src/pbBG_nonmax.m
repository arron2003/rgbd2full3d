function pball = pbBG_nonmax(im,radius,norient)
% function [pb,theta] = pbCGTG(im,radius,norient)
% 
% Compute probability of boundary using CG and TG.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% April 2003
%
% Edited by Derek Hoiem, Jan 2006: reduce radius size for large images,
% do not apply non-maxima suppression, 4 orientations


if nargin<2, radius=0.01; end
if nargin<3, norient=4; end

% beta from logistic fits (trainBG.m)
if radius==0.01,
  beta = [ -3.6944544e+00  1.0261318e+00 ];
  fstd = [  1.0000000e+00  3.7408935e-01 ];
  beta = beta ./ fstd;
else
  error(sprintf('no parameters for radius=%g\n',radius));
end

[imh, imw, unused] = size(im);
radius = norm([320 240]) / norm([imh imw]) * radius; % make radius equivalent to [240 320] case


% get gradients
[bg,gtheta] = detBG(im,radius,norient);

% compute oriented pb
[h,w,unused] = size(im);
pball = zeros(h,w,norient);
for i = 1:norient,
  b = bg(:,:,i); b = b(:);
  x = [ones(size(b)) b];
  pbi = 1 ./ (1 + (exp(-x*beta')));
  pball(:,:,i) = reshape(pbi,[h w]);
end

% mask out 1-pixel border where nonmax suppression fails
% pb(1,:) = 0;
% pb(end,:) = 0;
% pb(:,1) = 0;
% pb(:,end) = 0;
