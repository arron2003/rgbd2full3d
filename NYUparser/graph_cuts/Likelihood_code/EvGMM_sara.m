% EvGMM.m  -- as EvalGMM but with  set of pixels rather than sets of coords

% returns the cost = -2*log p(z|class) for
% probability at each pixel i, j in im
% assume assignment as max. prob. (hard assignment for labels!)

function [eImagek,assignment] = EvGMM(data, p, m, C)

% Initialise

% data                 = imk(pixels,:)';
[nDimensions, nData] = size(data);
eImagek               = zeros(1, nData);
nClasses             = length(p);

% Compute error e
% e = - 2 * log( p(class | data) )

for i = 1:nClasses
  if (p(i) > 0)  
    CInvRoot = chol(inv(C{i}));
    LogDetC  = log(det(C{i}));
    cData    = data - m{i} * ones(1, nData);
    ei = CInvRoot * cData;
    ei = ei .* ei;
    ei = sum(ei);
    % ei = ei ./ (nDimensions/3);
    ei = ei + LogDetC * ones(1, nData);
    ei = ei - ((2 * log(p(i))) * ones(1, nData));
    e(i, :) = ei;
  else
    e(i,:) = ones(1,nData)*Inf; 
  end
end;

% Evaluate energy using min approx across components
if (nClasses>1)
    [emin,assignment] = min(e);
else
    emin = e; 
    assignment = ones(size(emin));
end

% eImagek(pixels) = emin;
eImagek = emin;






