%---------------------------------------
% ComputeGMM.m
%
% Do clustering for general n-dimensional data
% data - size: dimension x data (e.g. dimension = 3)
% init: 1 for initialisation, otherwise 0
% p (init=0): mixture coefficient
% m (init=0): mean of the Gaussian
% C (init=0): convariance of the Gaussian
% sigma(init=0):        
%
% options: 
%
% GMMOptions.nComponents: number of GMM of the model
% GMMOptions.iterations: iterates for EM or kmeans method 
% GMMOptions.init: Initialisation Option: 0 - smart, 1 - random
% GMMOptions.method: 'kmeans' or 'em' - type of method to compute GMM                           				    
% GMMOptions.covarType: 'full' or 'spherical' - type of Gaussian with covariance matrix varying 
% GMMOptions.min_var: minimum variance in the covariance matrix (each dimension)
%
%-----------------------------------------

function [p, m, C] = ComputeGMM_cr( data, options, init, p, m, C, sigma);

% global variables
global Error im im_full; % for error messages

%-----------------------------------------
% Possible Parameter settings:
% 2 Gaussians with sigma 0.07
% data = zeros(3,200);
% data(:,1:100) = (randn(3,100)*0.07+repmat([0.3;0.3;0.3],1,100));
% data(:,101:200) = (randn(3,100)*0.07+repmat([0.7;0.7;0.7],1,100));
% data(:,201:500) = (randn(3,300)*0.02+repmat([0.1;0.9;0.1],1,300));
% index = find(data>1)
% data(index) = 1;
% index = find(data<0)
% data(index) = 0;
% init = 1;
% p = 0;
% m = 0;
% C = 0;
% sigma = 0;
% options.nComponents = 10; 
% options.iterations = 2;
% options.init = 0; 
% options.method = 'kmeans';                      				    
% options.covarType = 'spherical'; % 'full';
% options.min_var = (1/255)^2.0;
%
%-----------------------------------------

% Adaptions
n = options.nComponents;
useEM = strcmp(options.method, 'em');
full = strcmp(options.covarType, 'full');
nIterations = options.iterations;

% % get the sigma from C (TODO: this is incorrect since we should get it as input)
if (~full & ~init)
  count = length(C);
  sigma = zeros(1,count);
  for i = 1:count;
    hm1 = C{i};
    sigma(i) = sqrt(hm1(1,1));
  end 
end % if ~full

if (init) % initialise the GMM
  if (full)
    [p, m, C] = initgmm(data, n, options);
  else % if ~full
    [p, m, sigma] = initgmmiso(data, n, options);
  end
end


for iteration = 1:nIterations
    
    %
    % EM
    %
    if (useEM)
        if (full)
            [p, m, C] = gmmem(data, p, m, C);
        else % if ~full
            [p, m, sigma] = gmmisoem(data, p, m, sigma);
        end; % if ~full
  
    else % if ~em
    %
    % k-means
    %
        if (full)
            if (options.LazySnapping)
                [p, m, C] = gmmkm_LazySnapping(data, p, m, C);
            else
                [p, m, C] = gmmkm(data, p, m, C);
            end          

        else % if ~full
            [p, m, sigma] = gmmisokm(data, p, m, sigma);
        end; % if ~full
    end; % else if ~em

    if (length(p) <= 1) break; end;
    
end; % for iteration = 1:nIterations  

% check for Errors
if (Error > 0)
    return; 
end

if (length(p) < 1) 
   disp('Sorry but the Gaussians are degenerated - Start again');
   Error = 1;
   im = im_full;
   imshow(im_full);
   errordlg('Sorry, but the Gaussians were degenerated; restart with cut out','Error: Segmentation');
   return; 
end

if (~full)
  nDimensions = length(m{1});
  C = [];
  for i = 1:length(sigma);
    C{i} = sigma(i)^2 * eye(nDimensions, nDimensions);
  end; % for i = 1:nSigma
end; % if ~full

    %-------------------------------------
    % TEST
    %-------------------------------------

%     if (~full)
%         count = length(C);
%         for hi1 = 1:count
%             hm1 = C{hi1};
%             if (size(hm1,1)>3)
%                 disp('BE CAREFUL - ADAPT COV MATRIX');
%                 hm1(4,4) = 100; % hm1(4,4)*5;
%                 hm1(5,5) = 100; % hm1(5,5)*5;
%                 C{hi1} = hm1;
%             end
%         end
%     end
    
    %-------------------------------------

%-------------------------------------

% PlotColourSpace_vs2(data, data, p, m, C, [], m, C);
% PlotColourSpace_vs2(ones(3,1), ones(3,1), p, m, C, [], m, C);
% pause