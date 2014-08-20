% scrap for location prior on wall, floor and ceiling
load('../config/splits.mat');
%% collect cool stuff
xleftwall = [];
zwall = [];
hceil = [];
hfloor = [];

for i=trainNdxs'
  fn = dir(sprintf('../mat/%d_*.mat', i));
  fn = fn.name;
  fprintf('%s\n', fn);
  load(sprintf('../mat/%s', fn));
  for k=1:numel(model.objects)
    o = model.objects{k};
    if strcmp(model.objects{k}.model.type, 'layout')
      for p=1:numel(o.model.surfaces)
        plane = model.objects{k}.model.surfaces{p}.plane;
        if abs(plane(1))>(1-tolPara*2)
          v = - plane(4) / plane(1); % x wall
          xwall = [xwall, v];
        elseif abs(plane(2))>(1-tolPara*2)
          v = - plane(4) / plane(2); 
          if v>0
            % ceiling
            hceil = [hceil, v];
          else
            % floor
            hfloor = [hfloor, v];
          end
        elseif abs(plane(3))>(1-tolPara*2)
          v = - plane(4) / plane(3); % z wall
          zwall = [zwall, v];
        end
      end
    end
  end
end

%% estimate truncated normal distribution
layoutLoc = struct;

t = -xwall(xwall<0); save('tmp.txt', 't', '-ASCII'); unix('/usr/local/bin/R < truncNorm.R --vanilla'); para = load('tmp.txt');
layoutLoc.leftwall.mu = -para(1); layoutLoc.leftwall.sigma=para(2); 

t = xwall(xwall>0); save('tmp.txt', 't', '-ASCII'); unix('/usr/local/bin/R < truncNorm.R --vanilla'); para = load('tmp.txt');
layoutLoc.rightwall.mu = para(1); layoutLoc.rightwall.sigma=para(2);

t = -zwall; save('tmp.txt', 't', '-ASCII'); unix('/usr/local/bin/R < truncNorm.R --vanilla'); para = load('tmp.txt');
layoutLoc.frontwall.mu = -para(1); layoutLoc.frontwall.sigma=para(2);

t = hceil; save('tmp.txt', 't', '-ASCII'); unix('/usr/local/bin/R < truncNorm.R --vanilla'); para = load('tmp.txt');
layoutLoc.hceil.mu = para(1); layoutLoc.hceil.sigma = para(2);

t = -hfloor; save('tmp.txt', 't', '-ASCII'); unix('/usr/local/bin/R < truncNorm.R --vanilla'); para = load('tmp.txt');
layoutLoc.hfloor.mu = -para(1); layoutLoc.hfloor.sigma = para(2);

save('../config/layoutLocation.mat', 'layoutLoc');