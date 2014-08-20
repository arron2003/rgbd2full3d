function [UnaryPotentials,settings,assignmentF,assignmentB] = getLikelihood(original_image,settings,trimap,indexF,indexB)


options.nComponents = settings.GMM_nComponents; 
options.iterations = 4;
options.init = 0; %Initialisation Option: 0 - smart, 1 - random
options.method = 'kmeans';                      				    
options.covarType = 'full'; %'spherical';%
options.min_var = (1/255)^2.0;
options.LazySnapping = 0; 


% options.nComponents = 30; 
% options.iterations = 1;  % iterates the EM or kmeans method 
% options.init = 0; % 0 - smart, 1 - random
% options.method = 'kmeans';  % type of method to compute GMM                              				    
% options.covarType = 'full'; % type of Gaussian with covariance matrix varying 
% options.min_var = (1/255)^2.0; % minimum variance in the covariance matrix (each dimension)



%%

[M,N,S] = size(original_image);

imr = original_image(:,:,1);
img = original_image(:,:,2);
imb = original_image(:,:,3);

UnaryPotentials = zeros(M,N,2);

% Unknown data
dataU = zeros(3,size(trimap.U,1));
dataU(1,:) = imr(trimap.U);
dataU(2,:) = img(trimap.U);
dataU(3,:) = imb(trimap.U);

%Used to compute Foreground color model
dataF = zeros(3,size(indexF,1));
dataF(1,:) = imr(indexF);
dataF(2,:) = img(indexF);
dataF(3,:) = imb(indexF);

%Used to compute Background color model
dataB = zeros(3,size(indexB,1));
dataB(1,:) = imr(indexB);
dataB(2,:) = img(indexB);
dataB(3,:) = imb(indexB);

if (isfield(settings, 'pF'))
    [pF, mF, CF] = ComputeGMM_cr( dataF, options, 0, settings.pF, settings.mF, settings.CF, 0);
    [pB, mB, CB] = ComputeGMM_cr( dataB, options, 0, settings.pB, settings.mB, settings.CB, 0);
else
    [pF, mF, CF] = ComputeGMM_cr( dataF, options, 1, 0, 0, 0, 0);
    [pB, mB, CB] = ComputeGMM_cr( dataB, options, 1, 0, 0, 0, 0);
end

settings.pF = pF;
settings.mF = mF;
settings.CF = CF;
settings.pB = pB;
settings.mB = mB;
settings.CB = CB;

assignmentF = zeros(M,N);
assignmentB = zeros(M,N);

[fLikeVecU,assignmentF(trimap.U)] = EvGMM_sara(dataU, pF, mF, CF);
[bLikeVecU,assignmentB(trimap.U)] = EvGMM_sara(dataU, pB, mB, CB);

fLike = zeros(M,N);
bLike = zeros(M,N);

fLike(trimap.U) = fLikeVecU;
bLike(trimap.U) = bLikeVecU;



%% ------ Fixing labeled regions
k = 4*(settings.lambda1 + settings.lambda2)+10000;

if (settings.neighboorSystem)
    k = k + 4*(settings.lambda1 + settings.lambda2)/sqrt(2);
end
settings.k = k;


%k = 1000000;
%k=Inf;

bLike(trimap.B) = 0;
fLike(trimap.B) = k;
if (settings.fixed_foreground)
    %the foreground pixels are only fixed if the setting
    %fixed_foreground is one
    bLike(trimap.F) = k;
    fLike(trimap.F) = 0;
end

UnaryPotentials(:,:,1) = bLike;
UnaryPotentials(:,:,2) = fLike;

if(settings.border_around == 1)
    UnaryPotentials(1,:,2) = UnaryPotentials(1,:,2) + k;
    UnaryPotentials(end,:,2) = UnaryPotentials(end,:,2) + k;
    UnaryPotentials(2:end-1,1,2) = UnaryPotentials(2:end-1,1,2) + k;
    UnaryPotentials(2:end-1,end,2) = UnaryPotentials(2:end-1,end,2) + k;
end


%figure; imagesc(bLike-fLike);
%Likelihood_figure(gcf);

end