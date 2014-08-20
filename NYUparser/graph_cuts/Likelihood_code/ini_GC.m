
% Assumes variables 'original_image', 'user_input' and settings
% indexF and indexB are the pixels for training the color models
% structure trimap.F, trimap.B is used for hard_constraining fore and
% background
% trimap.U is the unknown data

[M,N,S] = size(original_image);

im_beta = Compute_beta(original_image,settings.neighboorSystem);

[UnaryPotentials,settings]  = getLikelihood(original_image,settings,trimap,indexF,indexB);

[LabelOut_GC, energy] = mex_graphCut(UnaryPotentials,settings.lambda1,settings.lambda2,settings.neighboorSystem,original_image,im_beta);

for k=1:settings.GrabCut_iterations - 1

    indexB = find(LabelOut_GC == 0);
    indexF = find(LabelOut_GC == 1);

    [UnaryPotentials,settings]  = getLikelihood(original_image,settings,trimap,indexF,indexB);
    
    [LabelOut_GC, energy_GC] = mex_graphCut(UnaryPotentials,settings.lambda1,settings.lambda2,settings.neighboorSystem,original_image,im_beta);
end

clear indexB indexF