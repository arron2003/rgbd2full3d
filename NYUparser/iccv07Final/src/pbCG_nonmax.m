function pball = pbCG_nonmax(im,radius,norient)
% function [pb,theta] = pbCGTG(im,radius,norient)
% 
% Compute probability of boundary using CG and TG.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% April 2003
%
% Edited by Derek Hoiem, Jan 2006: reduce radius size for large images,
% do not apply non-maxima suppression, 4 orientations


if nargin<2, radius=0.02; end
if nargin<3, norient=4; end % was 8

% beta from logistic fits (trainCG.m)
if all(radius==0.02),
  beta = [ -2.9216153e+00  2.1939403e-01  5.3764451e-01 ];
  fstd = [  1.0000000e+00  1.4210176e-01  1.9449891e-01 ];
  beta = beta ./ fstd;
else
  error(sprintf('no parameters for radius=%g\n',radius));
end
[imh, imw, unused] = size(im);
radius = norm([320 240]) / norm([imh imw]) * radius; % make radius equivalent to [240 320] case


% get gradients
[cg,theta] = detCG(im,radius,norient);
%[cg,tg,gtheta] = detCGTG(im,radius,norient);

% compute oriented pb
[h,w,unused] = size(im);
pball = zeros(h,w,norient);
for i = 1:norient,
  a = cg(:,:,2,i); a = a(:);
  b = cg(:,:,3,i); b = b(:);
  x = [ones(size(b)) a b];
  pbi = 1 ./ (1 + (exp(-x*beta')));
  pball(:,:,i) = reshape(pbi,[h w]);
end

% mask out 1-pixel border where nonmax suppression fails
% pb(1,:) = 0;
% pb(end,:) = 0;
% pb(:,1) = 0;
% pb(:,end) = 0;
