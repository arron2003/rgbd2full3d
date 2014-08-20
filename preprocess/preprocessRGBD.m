function data = processRoom(rgb, depth, rawDepths, opt)
  data.images = im2double(rgb); data.depths = depth; data.rawDepths = rawDepths;
  if opt.debug, fprintf('Preprocessing RGBD data..\n'); end
  p3d = rgb_plane2rgb_world_cropped(depth);
  p3d = reshape(p3d, [size(data.images,1), size(data.images,2), 3]);
  data.X = p3d(:,:,1);  data.Y = p3d(:,:,2);  data.Z = -p3d(:,:,3);
  
  if opt.debug, figure(1); imshow(data.images); drawnow;
    fprintf('Computing Normal..\n'); end
  [normalim, normalconf] = computeNormal(data);
  data.normal = normalim;
  
  if opt.debug, figure(2); imshow(abs(data.normal)); drawnow;
    fprintf('Computing Room Rotation..\n'); end
  data.R = computeRoomDirections(data);
  data.normal = reshape(...
  reshape(data.normal, [size(data.normal, 1)*size(data.normal, 2), 3])*data.R, ...
  [size(data.normal,1), size(data.normal,2), 3]);
  p3d = reshape([data.X(:) data.Y(:) data.Z(:)]*data.R, [size(data.normal,1), size(data.normal,2), 3]);
  data.X = p3d(:,:,1);  data.Y = p3d(:,:,2);  data.Z = p3d(:,:,3);
  
  if opt.debug, 
    r = [1 0 0; 0 1 0; 0 0 1]*data.R';
    figure(3), clf; hold on;
    plot_arrow(0, 0, r(1,1), r(1,2), 'linewidth', 3, 'color', [1 0 0], 'edgecolor', [1 0 0], 'headwidth', .05/norm(r(1,1:2)), 'headheight', .1/norm(r(1,1:2)));
    plot_arrow(0, 0, r(2,1), r(2,2), 'linewidth', 3, 'color', [0 1 0], 'edgecolor', [0 1 0], 'headwidth', .05/norm(r(2,1:2)), 'headheight', .1/norm(r(2,1:2)));
    plot_arrow(0, 0, r(3,1), r(3,2), 'linewidth', 3, 'color', [0 0 1], 'edgecolor', [0 0 1], 'headwidth', .05/norm(r(3,1:2)), 'headheight', .1/norm(r(3,1:2)));
    axis([-1 1 -1 1]); axis equal; axis off;
    set(gcf, 'Color', [0 0 0]); drawnow;
  end
  
end