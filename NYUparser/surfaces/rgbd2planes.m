% Finds major scene surfaces (planes) and rotates the entire scene so that
% the floor surface normal points in the same direction as the positive
% Z-axis in the right-handed coordinate frame (up).
%
% Args:
%   imgRgb - 480x640x3 uint8 image.
%   imgDepthOrig - 480x640 depth image (not cropped).
%   imgDepthRawOrig - 480x640 projected depth image, pre-In-painting.
%   imgNormals - HxWx3 image of surface normals.
%   normalConf - the confidence of the surface normal calculation at each
%                pixel location.
%
% Returns:
%   planeData - a struct containing information regarding the scene's
%               3d point cloud, surface normals and major scene surfaces.
function planeData = rgbd2planes(imgRgb, imgDepthOrig, imgDepthRawOrig, ...
    imgNormals, normalConf)
  
  assert(isa(imgRgb, 'uint8'));

  Consts;
  
  [mask, sz] = get_projection_mask();
  
  points3d = rgb_plane2rgb_world_cropped(imgDepthOrig);
  %points3d = points3d(mask, :);
  %keyboard;
  
  Xrgb = reshape(points3d(:,1), sz);
  Yrgb = reshape(points3d(:,2), sz);
  Zrgb = reshape(points3d(:,3), sz);
  
  % Need this for the maxDepth variable.
  camera_params;
  measuredZ = imgDepthRawOrig > 0 & imgDepthRawOrig < maxDepth;
  %measuredZ = reshape(measuredZ(mask), sz);

  normals = reshape(imgNormals, [prod(sz) 3]);

  % Find major planes using RANSAC
  [planes, pts_idx] = xyz2planes_ransac(Xrgb, Yrgb, Zrgb, normals, measuredZ);

  % Align to room coordinates
  Rf = get_room_directions(im2double(imgRgb), [Xrgb(:) Yrgb(:) Zrgb(:)], measuredZ, normals(measuredZ .* normalConf > 0.9, :));
  
  if all(Rf(:) == 0)
    Rf = eye(3);
  end
  
  % Make sure we don't flip the entire room...
  if Rf(1,1) < 0
    Rf = [-1 0 0; 0 1 0; 0 0 1] * Rf;
  end
  
  if Rf(2,2) < 0
    Rf = [1 0 0; 0 -1 0; 0 0 1] * Rf;
  end
  
  if Rf(3,3) < 0
    Rf = [1 0 0; 0 1 0; 0 0 -1] * Rf;
  end
  
  assert(Rf(1,1) >= 0);
  assert(Rf(2,2) >= 0);
  assert(Rf(3,3) >= 0);
  
%   K_rgb = struct('fx', fx_rgb, 'fy', fy_rgb, 'u0', cx_rgb, 'v0', cy_rgb);
  
  [Xf, Yf, Zf, normalsf] = rotate_scene(Rf, Xrgb, Yrgb, Zrgb, normals, planes);

  % Perform graph cut to reassign pixels to planes
  depthWeightInterpolated = 0.25;
  depthWeight = double(measuredZ);
  depthWeight(~measuredZ) = depthWeightInterpolated;

  
  pts3d = [Xf(:) Yf(:) Zf(:) ones(prod(sz), 1)];
  [planemap, planesf] = cut_planes_registered(im2double(imgRgb), ...
      pts_idx, pts3d, normalsf, Zrgb, depthWeight);

  pts3d = pts3d(:, 1:3);

  supportSurfaces = get_support_surfaces(pts3d, normalsf, sz(1), sz(2));

  % Infer the maximum extent of each plane and determine whether the plane is
  % a scene boundary, vertical, horizontal, etc.
  extrapolatedPlanesMask = extrapolatePlanesGC(planesf, planemap, pts3d);
  planeInfo = interpret_planes(planesf, extrapolatedPlanesMask);
  % drawExtrapolatedPlanes(planesf, extrapolatedPlanesMask, planeInfo, min(pts3d(:, 1:3)), max(pts3d(:, 1:3)), Pf, 10);
  % drawnow;
  % 
  % % Create an overhead view of the visible portion of the scene
%   settings = struct('gridsize', [500 500], 'scale', 5/500);
%   [floormaps, heightmap, gridX, gridZ] = getFloorMap(planesf, planeInfo, planemap, ...
%       extrapolatedPlanesMask, Pf, Xf, Yf, Zf, true(sz), settings);
  % 
  % figure(3), imagesc(heightmap), axis image, colormap jet
  % figure(4), imagesc(planemap), axis image, colormap jet
  %XYZim = cat(3, Xf, Yf, Zf);


  planeData = struct('planeMap', uint16(planemap), 'planeEqs', planes, ...
        'points3d', single(pts3d), 'normals', single(normals), 'R', Rf);    

  planeData.planeTypes = planeInfo;
  planeData.supportSurfaces = supportSurfaces;
end