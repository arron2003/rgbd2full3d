function [mode, score] = mean_shift(data, NUM_OF_TRAILS, h, epsilon, dist2)
% locating the mode of data (local maxima of density function)
% using mean shift.
% h - Gaussian kernel with parameter h (window size)
% eps - epsilon accuracy of mean shift step. i.e. how small m_hG to be
%       considered as no improvement - NOT RECCOMAND eps = 0!!
% data - a dxN matrix of N vector of d dimension
% dist2 - pointer ro function that measures distance square
%         the distance function should be able to take two d by N matrices
%         and return a vector of length N, where Ni is the distance^2
%         between the two i-th vectors.
%
%
% Copyright (c) Bagon Shai
% Department of Computer Science and Applied Mathmatics
% Wiezmann Institute of Science
% http://www.wisdom.weizmann.ac.il/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% 1. The above copyright notice and this permission notice shall be included in 
%     all copies or substantial portions of the Software.
% 2. No commercial use will be done with this software.
% 3. If used in an academic framework - a proper citation must be included.
%
% The Software is provided "as is", without warranty of any kind.
%
% Mar. 2007
%


%NUM_OF_TRAILS = 1;

[d N] = size(data);

NUM_OF_TRAILS = min(NUM_OF_TRAILS,N);
% choose the random starting points
rsp  =randperm(N);
rsp = rsp(1:NUM_OF_TRAILS);
%rsp = randsample(1:N,NUM_OF_TRAILS);
eps2 = epsilon*epsilon;


%x = data(:, rsp(1));
%[mode, score] = mean_shift_iteration(x, data, h, eps2, dist2);

mode = zeros(d, NUM_OF_TRAILS);
score = zeros(1, NUM_OF_TRAILS);
for ii = 1:NUM_OF_TRAILS
    x = data(:,rsp(ii));  % start from a random x
    [mode(:, ii), score(ii)] = mean_shift_iteration(x, data, h, eps2, dist2);
            
    %if ( s > score )
    %    mode = m;
    %    score = s;
    %end
end 

[sv, si] = sort(score);
score = sv;
mode = mode(:, si);
keep = true(NUM_OF_TRAILS, 1);
for ii = 1:NUM_OF_TRAILS-1
    d = dist2(repmat(mode(:, ii), [1 NUM_OF_TRAILS-ii]), mode(:, ii+1:end));
    if any(d < h^2)
        keep(ii) = false;
    end
end

score = score(keep);
mode = mode(:, keep);
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [mode, score] = mean_shift_iteration(x, data, h, eps2, dist2, per)
% iterate mean shift search from starting point x

% compute the first mean shift step
[ms, score] = m_hG(x, data, h, dist2);

% iterate untill no step is achieved
while  ( dist2(ms , x) >= eps2 ) 

    x = ms;
    [ms, score] = m_hG(x, data, h, dist2);
end

mode = ms;

% prune non local maxima points
% %     if ( ~ pertrube(x, data, h, eps2, dist2) )
% %         score = -1;
% %     end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LocalMax] = pertrube(x, data, h, eps2, dist2)
% we reahed a stationary point, we want to pertrube x a little
% and check if its a local maxima or not
% the critiria: move x a little if the mean shift returns to the same point
% its a local maxima.

% pertrubing by random vector of norm h/4
per = rand( size(data,1) ,1)-0.5;
per = h * per ./ ( sqrt( dist2(per, zeros(size(per)) ) ) );
[mode, score] = mean_shift_iteration(x + per, data, h, eps2, dist2);

% logical variable
% the mean shift step must be toward x.
LocalMax = ( dist2(mode, x) <= sqrt(eps2) );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ms, score]=m_hG(x, data, h, dist2)
% compute the mean shift of the data
% x - relative to x location
% data - dxN 
% h - gaussian window

%       SUM Xi exp( -.5 || dist(X,Xi)/h ||^2 )
% m   = ---------------------------------
%  h,G  SUM exp( -.5 || dist(X,Xi)/h ||^2 )

% e is the argument of the exponent in the term above: i.e. 1xN vector with
% ||X-Xi/h||^2 at the i-th column

e = dist2( repmat(x, 1, size(data,2)), data ) ./ (h.^2);
e = exp(-0.5*e);
score = sum(e); 

% the mean shift
ms =  (data*e' ./ score);

