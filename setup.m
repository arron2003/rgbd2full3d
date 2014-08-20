% Download basic data
if ~exist('processed_data')
  !wget -r -l1 -np "http://aqua.cs.uiuc.edu/processed_data" -P . -A "*.mat" --no-host-directories --no-verbose
end

if ~exist('mat')
  !wget -r -l1 -np "http://aqua.cs.uiuc.edu/mat" -P . -A "*.mat" --no-host-directories --no-verbose
end

!mkdir cache
!mkdir cache_gt
!mkdir cache_prop