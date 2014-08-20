% Estimate a rotation matrix so that the floor normal is in the Y direction
% and the wall normals are in the X and Z directions
%
% 1) Find straight lines and their XYZ directions
% 2) Find modes of surface normals, compute VP of normal and vanishing line
% of planes
% 3) Generate candidates from surface normal modes and straight lines
% 4) Score each candidate
%
% Input
%  imgRgb - HxWx3 matrix (double type)
%  points3d - Nx3 matrix of 3D points where N=HxW
%  imgValid - HxW logical matrix, true if 3d point at that pixel is valid
%  norms3d - Nx3 matrix of surface normals at each point.
function R = get_room_directions(imgRgb, points3d, imgValid, norms3d)

  straightnessThreshold = 30; % ratio of first to second eigenvalue
  minLen = 15; % minimum length of a line
  w_line = 0.3; %0.3;
  w_surf = 0.7; %0.7;

  imgGray = rgb2gray(imgRgb);

  % get straight line pixels and parameters
  [~, lineidx] = get_straight_line_segments(imgGray, minLen);

  % Start by getting a bunch of hypothesizes lines and planes in 3d.
  xyz_lines = get_lines_3d(lineidx, points3d, imgValid, straightnessThreshold);
  xyz_surf = get_surfaces_3d(norms3d);
  xyz_cand = [xyz_surf ; xyz_lines];
  
  ty = xyz_cand(abs(xyz_cand(:, 2))>0.8, :);
  Ny = zeros(size(ty, 1)*10, 3); 
  Nx = zeros(size(ty, 1)*10, 3); 
  Nz = zeros(size(ty,1)*10, 3); 

  maxpery = 10;
  scores = zeros(size(ty, 1)*maxpery, 1);
  n = 0;
  for k = 1:size(ty,1)
    distY = abs(xyz_cand*ty(k, :)');
    ind = find(distY<0.01);
    if numel(ind) > maxpery
      rp  =randperm(numel(ind));
      ind = ind(rp(1:maxpery));
    end
    for k2 = 1:numel(ind) 
      n=n+1;
      Ny(n, :) = ty(k, :);        
      v2 = cross(xyz_cand(ind(k2), :), Ny(n, :));
      v1 = cross(Ny(n, :), v2);
      if v1(1)^2>v2(1)^2
        Nx(n, :) = v1;
        Nz(n, :) = v2;
      else
        Nx(n, :) = v2;
        Nz(n, :) = v1;
      end

      distX = xyz_lines*Nx(n, :)';
      distY = xyz_lines*Ny(n, :)';
      distZ = xyz_lines*Nz(n, :)';        
      score_line = sum(exp(-distX.^2/0.01.^2)) + sum(exp(-distY.^2/0.01.^2)) + ...
          sum(exp(-distZ.^2/0.01.^2)); 

      distX = norms3d*Nx(n, :)';
      distY = norms3d*Ny(n, :)';
      distZ = norms3d*Nz(n, :)';        
      score_surf = sum(exp(-distX.^2/0.01.^2)) + sum(exp(-distY.^2/0.01.^2)) + ...
          sum(exp(-distZ.^2/0.01.^2));     

      scores(n) = score_line / size(xyz_lines, 1) * w_line + ...
          score_surf / size(norms3d, 1) * w_surf;
    end
  end

  % make it so floor points up
  if Ny(2)>0
      Ny = -Ny;
  end

  [mv, mi] = max(scores);
  R = [Nx(mi, :) ; Ny(mi, :) ; Nz(mi, :)];
  %disp(num2str(R));
end

function [xyz_lines, imgLines] = get_lines_3d(lineidx, points3d, imgValid, straightnessThreshold)
  [H, W] = size(imgValid);
  imgLines = zeros(H, W);

  numLines = numel(lineidx);

  xyz_lines = ones(numLines, 3);
  n = 0;
  for k = 1 : numLines

      ind = lineidx{k}(imgValid(lineidx{k}));
      if numel(ind) < 5
        continue;
      end

      XYZ = points3d(ind, :);
      XYZ = XYZ - repmat(mean(XYZ, 1), [size(XYZ, 1) 1]);
      %A = [XYZ ones(size(XYZ,1), 1)];
      [eigv, eigl] = eig(XYZ'*XYZ);    
      if eigl(3,3) / eigl(2,2) < straightnessThreshold
          continue;
      end
      n=n+1;
      N = eigv(1:3, 3); % / sqrt(sum(eigv(1:3, 1).^2));
      if max(abs(N))>max(N)
          N = -N;
      end

      %vp_lines(n, 1) = K.fx * N(1) / N(3) + K.u0;
      %vp_lines(n, 2) = K.fy * N(2) / N(3) + K.v0; 
      xyz_lines(n, :) = N / sqrt(sum(N.^2));
      if 0 
          if abs(xyz_lines(n, 2))>0.9
              imgLines(lineidx{k}) = 1;
          elseif abs(xyz_lines(n, 1)) >0.9
              imgLines(lineidx{k} + H * W) = 1;        
          elseif abs(xyz_lines(n, 3)) >0.9
              imgLines(lineidx{k} + 2 * H * W) = 1;        
          end
      end
      imgLines(ind) = n;
      % vp_lines(3) = 1 by initialization
  end

  xyz_lines = xyz_lines(1:n, :);
end

function xyz_surf = get_surfaces_3d(norms3d)

  ntrails = 100;
  h=0.01;
  epsilon = 0.001;
  thresh = 50;
  [mode, score] = mean_shift(norms3d', ntrails, h, epsilon, @meanshift_dist2);
  mode = mode(:, score>thresh)';
  xyz_surf = zeros(size(mode));
  for k = 1:size(mode, 1)
      [mv, mi] = max(abs(mode(k, :)*norms3d'));
      xyz_surf(k, :) = norms3d(mi, :);
  end
end
