function [bndinfo, pbim, gconf] = im2boundariesTopLevel(im, thresh, gcbase)
%
% [bndinfo, pbim, gconf] = im2boundariesTopLevel(im, thresh)
%
% High level function for estimating boundaries in an image using Hoiem et
% al. 2007 occlusion reasonining method.  
%
% Inputs:
%   im:       RGB double format image
%   thresh:   threshold values for segmentation hierarchies (default =
%             [0.105 0.25 0.6], threshold is lowest boundary confidence 
%             required not to merge to regions in hierarchical segmentation
%   gcbase:   used to set geometric context results directory and 
%             filenames if source code not available
%
% Outputs:
%   bndinfo.wseg: the superpixels
%   bndinfo.edges.indices: the indices for each edge
%   bndinfo.result.edgeProb: probability of edge being on
%   bndinfo.result.geomProb: probabilty of geometric label (gnd, vert, sky)
%   bndinfo.result.boundaries: most likely ege label (0=off, 1=on)
%   pbim: probability of boundary (Pb) image (imh, imw, 4 orient) 
%   gconf: surface likelihoods (imh, imw, [support, vert-planar/L/C/R,
%          non-planar solid/porous])
% Notes:
%   - length of result.boundaries is twice length of edges.indices.  If
%   label in first half is "on", left side occludes; if label in second
%   half is "on", right side occludes.   
%   - if geometric context source not available, no re-estimate of
%   geometric context based on 3rd iteration segmentation will occur. Note
%   that the benefits of this step to the occlusion reasoning are not
%   great.
%
% Citation:
%   D. Hoiem, A.N. Stein, A.A. Efros, M. Hebert, "Recovering Occlusion
%   Boundaries from a Single Image", ICCV 2007

if ~exist('thresh', 'var') || isempty(thresh)
    thresh = [0.105 0.25 0.6];
end

% load classifiers
load('boundaryClassifiers.mat');
load('continuityClassifiers.mat');
gclassifiers1 = load('ijcvClassifier.mat');
gclassifiers2 = load('perfectSegClassifierCv.mat');


% load geometric confidences
disp('geometry')
if exist('ijcvTestImage')
    imwrite(im, './tmpim.ppm');
    outfn = './tmpim.pnm';
    syscall = './segment 0.8 100 100 ./tmpim.ppm ./tmpim.pnm';
    [pg, tmp, imsegs] = ijcvTestImage(im, {syscall, outfn}, gclassifiers1);
    delete('./tmpim.ppm');
    delete('./tmpim.pnm');
    clear tmp
    gdata.imsegs = imsegs;
    gconf = pg2confidenceImages(imsegs, {pg});
    gconf = gconf{1}(:, :, 1:7);    
else
    if ~exist('gcbase', 'var')
        error('gcbase must be set if ijcv geometric context source code not available');
    end
    gcends = {'000', '090-045', '090-090', '090-135', '090-por', '090-sol', 'sky'};
    gconf = zeros([size(im, 1) size(im, 2) 7]);
    for k = 1:7        
        gconf(:, :, k) = im2double(imread([gcbase '.' gcends{k} '.pgm']));
    end
    gdata.imsegs = [];
    gclassifiers2 = [];
end

% create pb confidences
disp('pb')
pbim = pbCGTG_nonmax(im);


% get occlusion boundary labels
disp('boundaries')
bndinfo = im2boundaries(im, pbim, gconf, dtBnd, dtBnd_fast, dtCont, ...
    gdata, gclassifiers2, thresh);
