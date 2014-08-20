function [bndinfo, pball] = im2superpixels(im, initseg)

  if isa(im, 'uint8')
    im = im2double(im);
  end

  % NCS: running into numerical issues here.
  tmp = imfilter(im, fspecial('gaussian', 13, 2));
  tmp(tmp > 1) = 1;
  tmp(tmp < 0) = 0;

  pball = pbCGTG_nonmax(tmp);
  %pball(pball<0.025) = 0;
  wseg = pb2wseg(pball, 4000);
  wseg = double(wseg);

  imgEdges = wseg == 0;

  % force consistency wtih initseg
  if exist('initseg', 'var') && ~isempty(initseg)
      wseg = wseg + max(wseg(:))*(initseg-min(initseg(:)));
      [us, ui, uj] = unique(wseg(:));
      wseg = reshape(uj-1, size(im, 1), size(im, 2));
      wseg(imgEdges) = 0;
  end

  [edges, ~, neighbors, wseg] = seg2fragments(double(wseg), im, 25);
  bndinfo = processBoundaryInfo(wseg, edges, neighbors);
  
  % Just doing this to be consistent with the rest of the codebase. Better
  % to fix it here than in processBoundaryInfo that is part of Derek's
  % iccv07 code.
  bndinfo.imgRegions = wseg;
  bndinfo = rmfield(bndinfo, 'wseg');
end
