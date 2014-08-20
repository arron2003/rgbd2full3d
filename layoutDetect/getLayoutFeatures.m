function [label, feat, v] = getLayoutFeatures(plane, data, opt)
  %% extract layout feature for a plane
      [~, ass] = max(data.pm, [], 3);
      pfloor = data.pm(:,:,1).*double(ass==1);
      pwall = data.pm(:,:,2).*double(ass==2);
      pceiling = data.pm(:,:,3).*double(ass==3);
      pobject =  data.pm(:,:,4).*double(ass==4);
      
      if abs(plane(1))>(1-opt.tolPara*2)
        v = - plane(4) / plane(1); % facing x
        mask = normpdf((data.X-v)./data.Z, 0, opt.tolErr) .* normpdf(abs(data.normal(:,:,1)), 1, opt.tolPara) ...
          ./ normpdf(0, 0, opt.tolErr) ./ normpdf(0, 0, opt.tolPara);
        member = (data.X-v)./data.Z*sign(v);
        if v<0, s = 'leftwall'; else s = 'rightwall'; end
      elseif abs(plane(2))>(1-opt.tolPara*2)
        v = - plane(4) / plane(2); % facing y
        mask = normpdf((data.Y-v)./data.Z, 0, opt.tolErr) .* normpdf(abs(data.normal(:,:,2)), 1, opt.tolPara) ...
          ./ normpdf(0, 0, opt.tolErr) ./ normpdf(0, 0, opt.tolPara);
        member = (data.Y-v)./data.Z*sign(v);
        if v<0, s = 'hfloor'; else s = 'hceil'; end
      elseif abs(plane(3))>(1-opt.tolPara*2)
        v = - plane(4) / plane(3); % facing z
        mask = normpdf((data.Z-v)./data.Z, 0, opt.tolErr) .* normpdf(abs(data.normal(:,:,3)), 1, opt.tolPara) ...
          ./ normpdf(0, 0, opt.tolErr) ./ normpdf(0, 0, opt.tolPara);
        member = (data.Z-v)./data.Z*sign(v);
        s = 'frontwall';
      else
        label=[]; feat = []; v = []; return;
      end
      feat = zeros(12, 1);
      % features 1: 2D pixels count
      feat(1) = sum(mask(:));
      % features 2: wall pixel count
      feat(2) = sum(mask(:).*pwall(:));
      % features 3: floor structure pixel count
      feat(3) = sum(mask(:).*pfloor(:));
      % features 4: ceiling structure pixel count
      feat(4) = sum(mask(:).*pceiling(:));
      % features 5: object pixel count
      feat(5) = sum(mask(:).*pobject(:));
      % features 6: largest possible pixel count on this plane
      feat(11) = sum(member(:)<opt.tolErr);
      % normalized with respect to posible size
      feat(6:10) = (feat(1:5)+.1)/(10+feat(11));
      % features 12: location prior
      para = getfield(opt.layoutLoc, s);
      feat(12) = normpdf(v, para.mu, para.sigma) / ...
         ( normcdf(10, abs(para.mu), para.sigma) - normcdf(0, abs(para.mu), para.sigma) );
      label = s;
end