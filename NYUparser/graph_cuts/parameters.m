%parameters

settings.GMM_nComponents = 5; %5; %Number of GMM components
settings.GrabCut_iterations = 5; %5; %Number of times the color models are updated
settings.neighboorSystem = 1;
settings.p_norm = 1; % norm - note: it can also have multiple values. Example: settings.p = [1 2 3];
settings.width = 1; %Should be an odd number. it imposes a square width*width

settings.fixed_foreground = 1;
settings.border_around = 0;
settings.hard_constraints = 1;
settings.uniform_foreground = 0;


%image_file = 'examples\images\butterfly.jpg';
%trimap_file = 'examples\trimap\butterfly.bmp';
