function result = renderPovray(model)
  % flush the model into a obj file
  if isfield(model, 'vertices')
    fv = model;
  elseif isfield(model, 'R')
    data = model;
    coord = reshape(cat(3, data.X, data.Y, data.Z), [], 3);
    coord = coord * model.R';
    coord = reshape(coord, [size(data.X), 3]);
    fv=surf2patch(surf(coord(:,:,1), coord(:,:,2), coord(:,:,3)), 'triangles');
  elseif isfield(model, 'objects')
    fv.vertices=[];
    fv.faces=[];
    counter = 0;
    for k=1:numel(model.objects)
      for j=1:numel(model.objects{k}.mesh.comp)
        fv.vertices=[fv.vertices; model.objects{k}.mesh.comp{j}.vertices];
        fv.faces=[fv.faces; model.objects{k}.mesh.comp{j}.faces+counter];
        counter=counter+size(model.objects{k}.mesh.comp{j}.vertices,1);
      end
    end
    R = model.camera.R;
    fv.vertices = fv.vertices*R';
    fv.vertices(:,3) = fv.vertices(:,3);
  end
  
  fn = [fileparts( mfilename('fullpath') ) '/object.pov'];
  fid = fopen(fn, 'w');
  fwrite(fid, sprintf('mesh {\n'));
  s = '';
  for i=1:size(fv.faces, 1)
    s1 = sprintf('triangle { < %.5f, %.5f, %.5f >, < %.5f, %.5f, %.5f >, < %.5f, %.5f, %.5f > }\n', fv.vertices(fv.faces(i,:), :)');
    s = [s s1];
  end
  fwrite(fid, s);
  
  fwrite(fid, sprintf('}\n'));
  fclose(fid);
  
  fn = [fileparts( mfilename('fullpath') ) '/light.pov'];
  fid = fopen(fn, 'w');
  fwrite(fid, sprintf('#declare lightx=%f;\n#declare lighty=%f;\n#declare lightz=%f;\n', 0, max(0, max(fv.vertices(:,2))-.5), max(fv.vertices(:,3))/3));
  fclose(fid);
  
  %setenv('PATH', [getenv('PATH') ':/usr/local/bin:/home/guo29/local/bin']);
  ['/home/guo29/local/bin/povray +W561 +H427 +A0.2 ' fileparts( mfilename('fullpath') ) '/scene.pov']
  unix(['/home/guo29/local/bin/povray 2>/dev/null +W561 +H427 +A0.2 ' fileparts( mfilename('fullpath') ) '/scene.pov']);
  pause(0.2);
  result=im2double(imread([fileparts( mfilename('fullpath') ) '/scene.png']));
end
