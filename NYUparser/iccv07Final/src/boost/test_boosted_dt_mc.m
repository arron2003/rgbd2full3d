function confidences = test_boosted_dt_mc(classifier, features)
% confidences = test_boosted_dt_mc(classifier, features)
%
% Returns a log likelihod ratio for each class in the classifier    
% 
% Input:
%  classifier: boosted decision tree classifier
%  features:   classifier features (ndata, nvariables)
% Output:
%   confidences(ndata, nclasses): 
%      P(class=k|features) \propto 1./(1+exp(-confidences(k)))


if size(features, 2)~=classifier.wcs(1).dt(1).npred
    error('Incorrect number of attributes')
end

wcs = classifier.wcs;  
nclasses = size(wcs, 2);

ntrees = size(wcs, 1);

confidences = zeros(size(features, 1), nclasses);
for c = 1:nclasses    
    for t = 1:ntrees        
        if ~isempty(wcs(t,c).dt)            
            dt = wcs(t,c).dt;
            %nodes = treevalc(int32(dt.var), dt.cut, int32(dt.children(:, 1)), ...
            %        int32(dt.children(:, 2)), dt.catsplit(:, 1), features');  
            [class_indices, nodes, classes] = treeval(wcs(t, c).dt, features);        
            confidences(:, c) = confidences(:, c) + wcs(t, c).confidences(nodes);
        end        
    end
    confidences(:, c) = confidences(:, c) + classifier.h0(c);
end

   