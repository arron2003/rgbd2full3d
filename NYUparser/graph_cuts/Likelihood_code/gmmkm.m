% gmmkm.m

%
% Gaussian Mixture Model using K-means
%

function [p, m, C] = gmmkm(data, p0, m0, C0)

% global 
global Error im im_full; % for error messages

% Constants

[nDimensions, nData] = size(data);
nClasses             = length(p0);

minDetCov = max(1E-20, 10^(-200/nDimensions)); % minimum determinant of covariance matrix
% minDetCov = 0.0; % otherwise problems

if (nClasses <= 1)
  % Compute gaussian
  p(1) = 1;
  m{1} = mean(data')';
  C{1} = cov(data');
  return;
end;
  
% Compute error e
% e = - log( p(class | data) ) + const

for i = 1:nClasses
  [CInvRoot, err] = chol(inv(C0{i}));
  LogDetC  = log(det(C0{i}));
  cData    = data - m0{i} * ones(1, nData);
  ei = CInvRoot * cData;
  ei = 0.5 * sum(ei .^2);
  ei = ei + 0.5 * LogDetC;
  ei = ei - log(p0(i));
  e(i, :) = ei;
end;

% Assign class with mimimum error to each point
% and sort the classes

[minError, class]  = min(e);
[classSort, order] = sort(class);
classSortShiftR(2:nData + 1) = classSort(1:nData);
classSortShiftR(1) = 0;
classStart = find( classSortShiftR(1:nData) ~= classSort );

% Recompute class mean, covariance and prior

nClasses = length(classStart);
classStart(nClasses + 1) = nData + 1;

validClass = ones(1, nClasses);
for i = 1:nClasses
  classStarti   = classStart(i);
  classStartip1 = classStart(i + 1);
  nClass        = classStartip1 - classStarti;
  
  if (nClass < nDimensions) 
    % covariance degenerate
    validClass(1, i) = 0;
    continue;
  end;
  
  orderi        = order(classStarti : classStartip1 - 1);
  datai         = data(:, orderi);

  p(i)  = (classStartip1 - classStarti) / nData;
  m{i}  = mean( datai' )'; 
  C{i}  = cov( datai' );
  
  % test for valid covariance
  if (det(C{i}) < minDetCov)
    validClass(i) = 0;
  else
    [CInvRoot, err] = chol(inv(C{i}));
    if (err)
      validClass(i) = 0;
    end;
  end;

end;

validClassIndex = find(validClass);
if (size(validClassIndex,2) > 0)
    p = p(:, validClassIndex);
    p = p./(sum(p(:)));
    m = m(:, validClassIndex);
    C = C(:, validClassIndex);
else
   disp('No Gaussian Mixture Model could be determined - Start again');
   Error = 1;
   im = im_full;
   imshow(im_full);
   errordlg('Sorry, no Gaussian Mixture Model computable; restart with cut out','Error: Segmentation');
   p = zeros(1,1);
   m = zeros(1,1);
   C = zeros(1,1);
end


% TEST - show all assignments as a faked colour image
% if (nDimensions == 5)
% end
