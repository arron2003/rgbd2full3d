function [] = Interactive_connected_segmentation()


%%%%%%%%%%%%%%%%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_file = 'examples\images\35070.jpg';
trimap_file = 'examples\trimap\35070.bmp';

settings.lambda1 = 5;  %Ising Prior weight
settings.lambda2 = 95; %Edge sensitive weight
settings.GMM_nComponents = 5; %Number of GMM components
settings.GrabCut_iterations = 5; %Number of times the color models are updated
settings.neighboorSystem = 1;
settings.p_norm = 1; % norm - note: it can also have multiple values. Example: settings.p = [1 2 3];
settings.width = 1; %Should be an odd number. it imposes a square width*width

%%%%%%%%%%%%%%% End of Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting grab cut and updating color models');

addpath([cd '\Likelihood_Code']);

original_image = double(imread(image_file))/255;
user_input = imread(trimap_file);


settings.fixed_foreground = 1; 
settings.border_around = 0;
settings.hard_constraints = 1;
settings.uniform_foreground = 0;


[trimap] = ini_likelihood(user_input,settings );
indexF= trimap.F;
indexB = trimap.B;

ini_GC;

%Impose figure - width
w = ceil((settings.width-1)/2);
t = repmat(-w:1:w,[2*w+1,1]);
t_ = t';
f = [t(:) t_(:)];

disp('Computing Min-marginals');

%Compute min-Marginals
tic
[LabelOut_GC, energy_GC, min_marginals0_simple, min_marginals1_simple] =....
    mex_graphCut(UnaryPotentials,settings.lambda1,settings.lambda2,settings.neighboorSystem,original_image,im_beta,f);
disp(['Time to compute min-marginals - ' num2str(toc)]);

Plot_segmentation( original_image, LabelOut_GC);
title('Segmentation obtained using trimap. Click on the starting point and press enter');
[y, x] = getpts;

y= round(y);
x = round(x);

if size(x,1)>1
    warndlg('You selected more than one starting point. The first clique is taken as starting point');
end

y = y(1);
x = x(1);

if( y<1 || y > N || x < 1 || x > M)
    errordlg('The starting point is outside the image. Program stopped');
    return
end
    
start = x + (y-1) * M;

figure; imshow(original_image); title('Clique on the end points and press enter. You can choose multiple points.');
[y, x] = getpts;
y= round(y);
x = round(x);

diff_min_marginals = min_marginals1_simple - min(min_marginals1_simple(:));

for p = 1:length(settings.p_norm)
    p_norm = settings.p_norm(p);
    [parent_pointers,b,all_distances] = mex_shortest_path(diff_min_marginals.^p_norm,0,start-1,-1,0);

    paths_marginals = false(M,N);
    for k = 1:length(x)
        if( y(k)<1 || y(k) > N || x(k) < 1 || x(k) > M)
            continue;
        end
        final = x(k) + (y(k)-1) * M;
        path = zeros(size(parent_pointers));
        %path(final) = 1;
        path(get_figure(final,f,[M N])) = 1;
        parent = parent_pointers(final);
        while (parent > 0)
            path(get_figure(parent,f,[M N])) = 1;
            %path(parent) = 1;
            parent = parent_pointers(parent);
        end
        paths_marginals = paths_marginals | (path==1);
    end

    imr = original_image(:,:,1);
    img = original_image(:,:,2);
    imb = original_image(:,:,3);

    imr(paths_marginals) = 0;
    img(paths_marginals) = 0;
    imb(paths_marginals) = 0;

    imr(paths_marginals) = 1;
    figure; imshow(cat(3,imr,img,imb)); title(['Paths ' num2str(p_norm) '-norm']);

    a = UnaryPotentials(:,:,1);
    a(paths_marginals) = 100000;
    aux = UnaryPotentials;
    aux(:,:,1) = a;
    [LabelOut_GC, energy_GC] = mex_graphCut(aux,settings.lambda1,settings.lambda2,settings.neighboorSystem,original_image,im_beta);
    Plot_segmentation(original_image, LabelOut_GC); title('Segmentation');
end

end

function [f] = get_figure(pixel,fig,image_size)

    [x,y] = ind2sub(image_size,pixel);
    
    a = repmat([x y], [size(fig,1) 1]) + fig;
    
    ind = (a(:,1)>0 & a(:,1)<= image_size(1) & a(:,2)>0 & a(:,2)<= image_size(2));
    
    a = a(ind,:);
    
    f = sub2ind(image_size,a(:,1),a(:,2));

end