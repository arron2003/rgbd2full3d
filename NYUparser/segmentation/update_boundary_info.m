function [boundaryInfo2, err] = update_boundary_info(boundaryInfo, result, im)

  if isstruct(result) && isfield(result, 'regions')
    regions = result.regions;


    empty = false(numel(regions), 1);
    splab = zeros(boundaryInfo.nseg, 1);
    for k = 1:numel(regions)
        splab(regions{k}) = k;
        empty(k) = isempty(regions{k});        
    end
    regions(empty) = [];
  else
    splab = result;
    regions = cell(max(splab), 1);
    empty = false(numel(regions), 1);
    for k = 1:numel(regions)
        regions{k} = find(splab==k);
        empty(k) = isempty(regions{k});  
    end
    regions(empty) = [];
  end

  imgRegions = splab(boundaryInfo.imgRegions);

  % XXX only reading image until seg2fragments gets fixed
  %im = imread(['/usr1/projects/dhoiem/iccv07/iccvGroundTruth/images/' boundaryInfo.imname]);
  %im = im2double(im);

  [edges, juncts, neighbors, imgRegions] = seg2fragments(imgRegions, im, 1);
  boundaryInfo2 = processBoundaryInfo(imgRegions, edges, neighbors);

  if isfield(boundaryInfo, 'imname')
      boundaryInfo2.imname = boundaryInfo.imname;
  end

  %boundaryInfo2.imname = boundaryInfo.imname;

  stats = regionprops(boundaryInfo.imgRegions, 'Area');
  area = cat(1, stats(:).Area);
  boundaryInfo2.spArea = area;

  if isfield(boundaryInfo, 'type')
      boundaryInfo2.type = boundaryInfo.type;
      boundaryInfo2.names = boundaryInfo.names;
      [boundaryInfo2.labels, err] = iccvTransferLabels(boundaryInfo.labels, regions, area);
      boundaryInfo2 = processGtBoundaryLabels(boundaryInfo2);

  else
      err = nan;
  end

  boundaryInfo2 = orderfields(boundaryInfo2);
  boundaryInfo2.imgRegions = boundaryInfo2.wseg;
  boundaryInfo2 = rmfield(boundaryInfo2, 'wseg');
end