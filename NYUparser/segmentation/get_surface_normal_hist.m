function [incHists azHists] = get_surface_normal_hist(imgRegions, ...
    inclinations, binsInc, azimuths, binsAz)
  imgRegions = int32(imgRegions);
  inclinations = double(inclinations);
  binsInc = double(binsInc);
  azimuths = double(azimuths);
  binsAz = double(binsAz);
  
  [incHists, azHists] = mex_surface_normal_hist(imgRegions, ...
      max(imgRegions(:)), inclinations, binsInc, azimuths, binsAz);
    
  incHists = incHists';
  azHists = azHists';
end