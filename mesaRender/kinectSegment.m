function [I, D] = kinectSegment(fvc)
  
  fvc_org = fvc; 
  fvc2 = struct('faces', [], 'vertices', [], 'color', []);
  for i=1:numel(fvc_org)
    fvc1=fvc_org(i);
    fvc3.faces = reshape(1:size(fvc1.faces,1)*3, [size(fvc1.faces,2) size(fvc1.faces,1)])';
    t = reshape(fvc1.faces', [], 1);
    fvc3.vertices = fvc1.vertices(t, :);
    fvc3.color = ones(size(fvc3.vertices))*i/numel(fvc_org);
    
    fvc2.faces = cat(1, fvc2.faces, size(fvc2.vertices,1)+fvc3.faces);
    fvc2.color = cat(1, fvc2.color, fvc3.color);
    fvc2.vertices = cat(1, fvc2.vertices, fvc3.vertices);
  end
  fvc=fvc2;
  
  fvc.modelviewmatrix=create_lookat_matrix([-0.0244, -0.0259, 0], [-0.0244, -0.0259, 1], [0, 1, 0]);
  p = 7.8e-6;
  n = 519.46961112127485*p; f = 5;
  r = 1.0012*280*p; l = -1.0012*280*p; t = 213*p; b = -213*p;
  
  fvc.projectionmatrix=[2 * n / (r - l) 0 0 0; 0 2 * n / (t - b) 0 0; 0, 0, -(f + n) / (f - n), -(2 * f * n) / (f - n); 0 0 -1 0];
  fvc.vertices = double(fvc.vertices);
  fvc.viewport=[0 0 561 427];
  fvc.culling=0;
  fvc.lightposition=[0 0 0 0];
  fvc.enabledepthtest=1;
  I = zeros (561,427, 6);
  I(:,:,5)= inf;
  I=renderpatch(I,fvc);
  d = I(:,:,5);
  d = 2*n*f./(f+n-d*(f-n));
  D = flipud(d');
  I = flipdim(permute(I(:,:,1:3), [2 1 3]),1);
  I = int16(round(I(:,:,1)*numel(fvc_org)));
end
