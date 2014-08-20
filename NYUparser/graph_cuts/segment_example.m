
%%%%%%%%%%%%%%%%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imName = 'GT22';
scaling = 5; sF = 1/scaling;

image_file = sprintf('.\\matting_data\\inputData\\%s.bmp',imName);
trimap_file = sprintf('.\\matting_data\\scribbles\\%s.bmp',imName);

settings.lambda1 = 5;  %Ising Prior weight
settings.lambda2 = 90; %Edge sensitive weight

settings.GMM_nComponents = 5; %5; %Number of GMM components
settings.GrabCut_iterations = 5; %5; %Number of times the color models are updated
settings.neighboorSystem = 1;
settings.p_norm = 1; % norm - note: it can also have multiple values. Example: settings.p = [1 2 3];
settings.width = 1; %Should be an odd number. it imposes a square width*width

addpath([cd '\Likelihood_Code']);
original_image = imread(image_file);
user_input = imread(trimap_file);
%keyboard;
[h w c] = size(original_image);
nh = ceil(sF*h)* (1/sF);
nw = ceil(sF*w)* (1/sF);

original_image = double(imresize(original_image, [nh nw]))/255;
user_input = imresize(user_input, [nh nw]);

settings.fixed_foreground = 1;
settings.border_around = 0;
settings.hard_constraints = 1;
settings.uniform_foreground = 0;

[trimap] = ini_likelihood(user_input,settings );
indexF= trimap.F;
indexB = trimap.B;

[M,N,S] = size(original_image);
im_beta = Compute_beta(original_image,settings.neighboorSystem);
[UnaryPotentials,settings]  = getLikelihood(original_image,settings,trimap,indexF,indexB);

[gcseg, gc_energy] = mex_graphCut(UnaryPotentials,settings.lambda1,settings.lambda2,settings.neighboorSystem,original_image,im_beta);

subplot(3,1,1);
imagesc(original_image);
subplot(3,1,2);
imagesc(user_input);
subplot(3,1,3);
imagesc(gcseg);
