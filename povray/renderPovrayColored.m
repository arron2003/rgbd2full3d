function result = renderPovrayColored(data, opt)
  % flush the model into a obj file
  fn = [fileparts( mfilename('fullpath') ) '/object_colored.pov'];
  S = kinectSegment(data.final_fvc);
  im = reshape(im2double(data.images) ,[], 3);
  fvc = struct('vertices', [], 'faces', [], 'color', []);
  for ii=1:numel(data.final_fvc)
    c = mean(im(S==ii,:));
    if isnan(c(1)), c=[0 0 0]; end
    fvc(ii) = struct('vertices', data.final_fvc(ii).vertices, ...
      'faces', data.final_fvc(ii).faces, 'color', c);
  end
  
  fid = fopen(fn, 'w');
  for j=1:numel(fvc)
    s = '';
    fwrite(fid, sprintf('union{\nmesh {\n'));
    fv = fvc(j);
    for i=1:size(fv.faces, 1)
      s1 = sprintf('triangle { < %.5f, %.5f, %.5f >, < %.5f, %.5f, %.5f >, < %.5f, %.5f, %.5f > }\n', fv.vertices(fv.faces(i,:), :)');
      s = [s s1];
    end
    fwrite(fid, s);
    fwrite(fid, sprintf('}\n'));
    fwrite(fid, sprintf('pigment {\ncolor rgb<%f,%f,%f>\n}\nfinish \n{ ambient 0.3\n diffuse 0.5\n roughness 0.1\n}\n}\n', fv.color.^(1.2)));
  end
  fclose(fid);
  
  fn = [fileparts( mfilename('fullpath') ) '/light.pov'];
  fid = fopen(fn, 'w');
  fv .vertices= cat(1, fvc.vertices);
  fwrite(fid, sprintf('#declare lightx=%f;\n#declare lighty=%f;\n#declare lightz=%f;\n', 0, max(0, max(fv.vertices(:,2))-.5), max(fv.vertices(:,3))/3));
  fclose(fid);
  
  %setenv('PATH', [getenv('PATH') ':/usr/local/bin:/home/guo29/local/bin']);
  [opt.povray_bin ' +W561 +H427 +A0.2 ' fileparts( mfilename('fullpath') ) '/scene_colored.pov']
  cd('povray');
  unix([opt.povray_bin ' 2>/dev/null +W561 +H427 +A0.2 ' fileparts( mfilename('fullpath') ) '/scene_colored.pov']);
  cd('..');
  pause(0.2);
  result=im2double(imread([fileparts( mfilename('fullpath') ) '/scene_colored.png']));
  if opt.debug
    imshow(result);
  end
end
