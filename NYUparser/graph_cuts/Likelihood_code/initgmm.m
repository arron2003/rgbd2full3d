% initgmm.m

% Initialise mixture components evenly distributed along 
% 1st principal direction, take the middle value as mean

function [p, m, C] = initgmm(data, nClasses, options)

%---------------------------------------
% Input Settings

% save '..\Temp\initgmm.mat' data nClasses option;
% load '..\Temp\initgmm.mat';

%---------------------------------------

if (nClasses <= 1)
  % Compute gaussian
  p{1} = 1;
  m{1} = mean(data')';
  C1 = cov(data');
  
  % check cov.
  if (C1(1,1) < options.min_var )
    C1(1,1) = options.min_var;
  end
  if (C1(2,2) < options.min_var )
    C1(2,2) = options.min_var;
  end
  if (C1(3,3) < options.min_var )
    C1(3,3) = options.min_var;
  end
  
  C{1} = C1;
  
  return;
end;

[nDimensions, nData] = size(data);

CData = cov(data');
[V, D] = eig(CData); % eigenvalues D(iagonal matrix) and eigenvectors V (in column) as X*V = V*D. where D93,3) is the largest eigenvalue

v1 = V(:, nDimensions); % 1st principle component 
d1 = D(nDimensions, nDimensions); 

step         = nData / nClasses;
stepIndexTmp = step / 2;
for i = 1:nClasses
  stepIndex(i) = ceil(stepIndexTmp); % these are the middle number of each class n.
  stepIndexTmp = stepIndexTmp + step;
end;

[sortedData, sortIndex] = sort(v1' * data); % map to 1st principle component and sort it

for i = 1:nClasses
  p(i) = 1 / nClasses;
  if (options.init == 0)
     m{i} = data(:, sortIndex(stepIndex(i))); % this is the middle data value of the nth classe
  else
     m{i} = rand(nDimensions,1);
  end
  C1 = CData / nClasses^2;
  
  % check minimum cov
  if (C1(1,1) < options.min_var )
    C1(1,1) = options.min_var;
  end
  if (C1(2,2) < options.min_var )
    C1(2,2) = options.min_var;
  end
  if (C1(3,3) < options.min_var )
    C1(3,3) = options.min_var;
  end
  
  C{i} = C1;
  
end;
  

