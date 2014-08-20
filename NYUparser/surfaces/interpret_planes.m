% Labels planes as being vertical, horizontal, and/or boundary planes.
%
% Args:
%   planeEqs - Px4 matrix of plane equations where P is the number of total
%              planes.
%   planeMasks - HxWxP logical matrix of masks for each of the planes.
%
%
% Returns:
%   planeInfo - a struct containing indicator variables for the different
%               plane types.
function planeInfo = interpret_planes(planeEqs, planeMasks)
  assert(size(planeEqs,2) == 4);
  assert(isa(planeMasks, 'logical'));

  [H, ~, P] = size(planeMasks);

  isVertical = planeEqs(:, 2).^2 < 0.01;
  isHorizontal = planeEqs(:, 2).^2 > 0.9;
  isBoundary = false(P, 1);
  boundaryExtent = zeros(P, 2);

  % Check which planes are boundary planes: should be vertical or horizontal,
  % and at the borders of image should be either visible or occluded
  for pp = 1 : P    
    columns = sum(planeMasks(:,:,pp), 1) > 0.9 * H;
    if (isVertical(pp) || isHorizontal(pp)) && sum(columns)>0.1            
      isBoundary(pp) = true;
      extent = mean(planeMasks(:, :, pp), 1)>0.9;
      ind = find(extent);
      boundaryExtent(pp, 1) = ind(1);
      boundaryExtent(pp, 2) = ind(end);
    end    
  end

%   if 0 
%   % Adjust for overlapping boundary planes of similar orientations
%     for p=1:P   
%         if isBoundary(p)
%             ind = abs(planeEqs(:, 1:3)*planeEqs(p, 1:3)')>0.9;
%             ind(p) = false;
%             infront = false(P, 1);
%             infront(ind) = isBoundary(ind) & (planeEqs(ind, 4).^2<planeEqs(p, 4).^2);        
%             frontLeft = infront & (boundaryExtent(:, 2)<boundaryExtent(p, 2));
%             if any(frontLeft)
%                 boundaryExtent(p, 1) = max(boundaryExtent(p, 1), max(boundaryExtent(frontLeft, 2)));        
%             end        
%             frontRight = infront & (boundaryExtent(:, 1)>boundaryExtent(p, 1));
%             if any(frontRight)
%                 boundaryExtent(p, 2) = min(boundaryExtent(p, 2), min(boundaryExtent(frontRight, 1)));
%             end
%         end
%     end
%   end 

  planeInfo = struct('isVertical', isVertical, ...
                     'isHorizontal', isHorizontal, ...
                     'isBoundary', isBoundary, ...
                     'boundaryExtent', boundaryExtent);
end