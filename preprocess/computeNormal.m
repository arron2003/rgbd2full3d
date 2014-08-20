function [normalim, confim] = computeNormal(data, normalim)
flag_caliberation = false;
if nargin==2
  flag_caliberation = true;
end

bw = 1:15; % blockwidth
depththresh = 0.05; % relative depth threshold
depththresh2 = 0.0075; % a more accurate one, for caliberation
cradius = 10;
X = data.X; Y = data.Y; Z = data.Z;
depths = data.depths;

[imh, imw] = size(Z);
npix = imh*imw;

pts = [X(:) Y(:) Z(:) ones(npix, 1)];

[u, v] = meshgrid(1:imw, 1:imh);

N=8;
for i=1:N
  nu(i,:) = round((rand(1,2*N)-.5)*2*i^2);
  nv(i,:) = round((rand(1,2*N)-.5)*2*i^2);
end

% [gx, gy]=gradient(data.depths);
% s = sqrt((gx.^2+gy.^2))./(data.depths.^2)/.0075;

[nu, nv] = meshgrid([-bw 0 bw], [-bw 0 bw]);

nx = zeros(imh, imw);
ny = zeros(imh, imw);
nz = zeros(imh, imw);
nd = zeros(imh, imw);
confim = zeros(imh, imw);

if flag_caliberation
  p3d=cat(3, X, Y, Z);
  validimage = false([size(X, 1) size(X, 2) numel(nu(:))]);
  for j=1:numel(nu(:))
    s=(p3d-shift2d(p3d,[nu(j) nv(j)]));
    validimage(:,:,j) = abs(sum(normalim.*s, 3)) < depths.^2*depththresh2;
  end
  validimage = reshape(validimage, [size(X, 1)*size(X, 2) numel(nu(:))]);
end

ind_all = find(Z);
I = reshape(data.images, [size(data.images,1)*size(data.images,2) 3]);
for k = ind_all(:)'
    if flag_caliberation
      u2 = u(k)+nu(validimage(k,:));
      v2 = v(k)+nv(validimage(k,:));
    else
      u2 = u(k)+nu;
      v2 = v(k)+nv;
    end
    
    % check that u2 and v2 are in image
    valid = (u2 > 0) & (v2 > 0) & (u2 <= imw) & (v2 <= imh);
    u2 = u2(valid);
    v2 = v2(valid);
    ind2 = v2 + (u2-1)*imh;
    
    % check that depth difference is not too large
    %valid = min(abs(pts(ind2, 1:3)-pts(k, 1:3)), 2) < Z(k)*depththresh; %Z(ind2)-Z(k)) < Z(k)*depththresh;
    
    %valid = abs(depths(ind2)-depths(k)) < depths(k).^2*depththresh2;
    %valid = abs(depths(ind2)-depths(k)) < depths(k)*depththresh;
    valid = sqrt((X(ind2)-X(k)).^2+(Y(ind2)-Y(k)).^2+(Z(ind2)-Z(k)).^2) < depths(k).^2*depththresh2;
    valid2 = abs(I(ind2, 1) - I(k, 1))<cradius & abs(I(ind2, 2) - I(k, 2))<cradius & abs(I(ind2, 3) - I(k, 3))<cradius;
    valid = valid & valid2;
    u2 = u2(valid);
    v2 = v2(valid);

    ind2 = v2 + (u2-1)*imh;
    
    if numel(u2)<3, continue; end
    
    A = pts(ind2, :);        
    [eigv, l] = eig(A'*A);
    nx(k) = eigv(1,1);
    ny(k) = eigv(2,1);
    nz(k) = eigv(3,1);
    nd(k) = eigv(4,1);
    confim(k) = 1- sqrt(l(1) / l(2,2)); 
    
end

% normalize so that first three coordinates form a unit normal vector and
% the largest normal component is positive
planeim = cat(3, nx, ny, nz, nd); % ./ repmat(len, [1 1 4]);
len = sqrt(nx.^2+ny.^2+nz.^2);
planeim = planeim ./ repmat(len+eps, [1 1 4]);
normalim = planeim(:, :, 1:3);
normalim = -normalim.*repmat(sign(sum(cat(3, X, Y, Z).*normalim, 3)), [1 1 3]);

