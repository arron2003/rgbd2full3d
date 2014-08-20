function result = findICP(data, opt, bbox1, bbox2, coord1, reg1mask, org_fvc, canScale, canRotate, previous_result)
  result = struct('filename', {}, 'gtid', [], 'label', {},...
    'rotate', [], 'scaling', [], 'score', [], 'faces', [], ...
    'vertices', [], 'depthErr', []);
  depthError = inf;
  if canRotate, listRotation = pi/4*[0 1 -1 2 -2];
  else listRotation = 0; end
  if canScale, listScaling = 2.^[1/3 -1/3 2/3 -2/3 1 -1];
  else listScaling = 1; end

  if nargin>9
    result(1) = previous_result;
    depthError = previous_result.depthErr;
  end
  fvc = org_fvc;
  for s = listScaling
    for rotY = listRotation
      % inital object rendering
      v = bsxfun(@minus, org_fvc.vertices, bbox2.centroid);
      v = bsxfun(@times, v, s); 
      Ry=[cos(rotY) 0 sin(rotY);
          0 1 0;
          -sin(rotY) 0 cos(rotY)];
      v = v*Ry;
      v = bsxfun(@plus, v, bbox1.centroid); 
      % weird historical issue: change z , apply R and change z
      v = v*data.R'; fvc.vertices=v;
      
      % out of scope object
      if max(fvc.vertices(:,3))>0, continue; end
      D = kinectMesa(fvc);
      
      % object whose depthError is bounded
      if (abs(sum(~~D(:))-sum(reg1mask(:)))*opt.missing_C)>depthError
        continue;
      end

      % do ICP
      coord=rgb_plane2rgb_world_cropped(D); 
      coord(:,3) = -coord(:,3); coord = coord*data.R; 
      
      newp = coord( D~=0, :);
      coord1 = get_pcd_sample(coord1, 300, .9);
      try %possible failure if rendered point outside of the image
        newp = get_pcd_sample(newp, 300, .9);
      catch
        continue;
      end
      %[~, M1]=ICP_translation(newp, coord1, struct('Registration','Translational', 'Optimizer', 'closed', 'Verbose', 1));
      M = simpleICP(coord1, newp);
      newp = bsxfun(@plus, newp, M);
      %figure(2), clf; hold on; scatter3(newp(:,1), newp(:,2), newp(:,3), 'b'); scatter3(coord1(:,1), coord1(:,2), coord1(:,3), 'r');
      %keyboard;
      v = bsxfun(@plus, v, M);
      fvc.vertices = v;
      %[I, D] = kinectCamera(fvc);
      D = kinectMesa(fvc);
      
      m = ~~D;
      t = sum(abs(data.depths(m&reg1mask)-D(m&reg1mask))) ...
        + sum(max(data.depths(m&~reg1mask)-D(m&~reg1mask),0)) ...
        + opt.missing_C * sum(sum(~m&reg1mask));

      if t<depthError
        depthError = t;
        result(1).depthErr = t;
        result(1).rotate = rotY;
        result(1).scaling = s;
        result(1).vertices = fvc.vertices;
        result(1).faces = fvc.faces;
      end
    end
  end
end