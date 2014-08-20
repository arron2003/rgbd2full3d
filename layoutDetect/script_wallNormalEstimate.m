% get variance of wall-normal 
fnlist = dir('../processed_data/*.mat');
fnlist = {fnlist.name};
totPt = 0;
totDiv = 0;
for i = 1:numel(fnlist)
  load(['../processed_data/' fnlist{i}]);
  iswall = ismember(data.label_names, 'wall');
  iswall = [false iswall];
  mask = iswall(data.gt+1) & (data.rawDepths<3) & (data.rawDepths>.2);
  normalmax = max(abs(data.normal), [], 3);
  totDiv = totDiv + sum((normalmax(mask)-1).^2);
  totPt = totPt + sum(mask(:));
  fprintf('%s %.5f\n', fnlist{i}, sqrt(totDiv/totPt));
end
