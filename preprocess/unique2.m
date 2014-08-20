function varargout = unique_new(varargin) 
%UNIQUE Set unique.
%   C = UNIQUE(A) for the array A returns the same values as in A but with 
%   no repetitions. C will be sorted.    
%  
%   C = UNIQUE(A,'rows') for the matrix A returns the unique rows of A.
%   The rows of the matrix C will be in sorted order.
%  
%   [C,IA,IC] = UNIQUE(A) also returns index vectors IA and IC such that
%   C = A(IA) and A = C(IC).  
%  
%   [C,IA,IC] = UNIQUE(A,'rows') also returns index vectors IA and IC such
%   that C = A(IA,:) and A = C(IC,:). 
%  
%   [C,IA,IC] = UNIQUE(A,OCCURRENCE) and
%   [C,IA,IC] = UNIQUE(A,'rows',OCCURRENCE) specify which index is returned
%   in IA in the case of repeated values (or rows) in A. The default value
%   is OCCURENCE = 'first', which returns the index of the first occurrence 
%   of each repeated value (or row) in A, while OCCURRENCE = 'last' returns
%   the index of the last occurrence of each repeated value (or row) in A.
%  
%   [C,IA,IC] = UNIQUE(A,'stable') returns the values of C in the same order
%   that they appear in A, while [C,IA,IC] = UNIQUE(A,'sorted') returns the
%   values of C in sorted order. If A is a row vector, then C will be a row
%   vector as well, otherwise C will be a column vector. IA and IC are
%   column vectors. If there are repeated values in A, then IA returns the
%   index of the first occurrence of each repeated value.
% 
%   [C,IA,IC] = UNIQUE(A,'rows','stable') returns the rows of C in the same
%   order that they appear in A, while [C,IA,IC] = UNIQUE(A,'rows','sorted')
%   returns the rows of C in sorted order.
% 
%   The behavior of UNIQUE has changed.  This includes:
%     -	occurrence of indices in IA and IC switched from last to first
%     -	IA and IC will always be column index vectors
% 
%   If this change in behavior has adversely affected your code, you may 
%   preserve the previous behavior with:
% 
%        [C,IA,IC] = UNIQUE(A,'legacy')
%        [C,IA,IC] = UNIQUE(A,'rows','legacy') 
%        [C,IA,IC] = UNIQUE(A,OCCURRENCE,'legacy')
%        [C,IA,IC] = UNIQUE(A,'rows',OCCURRENCE,'legacy')
%
%   Examples:
%
%       a = [9 9 9 9 9 9 8 8 8 8 7 7 7 6 6 6 5 5 4 2 1]
%
%       [c1,ia1,ic1] = unique(a)
%       % returns
%       c1 = [1 2 4 5 6 7 8 9]
%       ia1 = [21 20 19 18 16 13 10 6]
%       ic1 = [8 8 8 8 8 8 7 7 7 7 6 6 6 5 5 5 4 4 3 2 1]
%
%       [c2,ia2,ic2] = unique(a,'stable')
%       % returns
%       c2 = [9 8 7 6 5 4 2 1]
%       ia2 = [1 7 11 14 17 19 20 21]'
%       ic2 = [1 1 1 1 1 1 2 2 2 2 3 3 3 4 4 4 5 5 6 7 8]'
%
%       c = unique([1 NaN NaN 2])
%       % NaNs compare as not equal, so this returns
%       c = [1 2 NaN NaN]
%
%   Class support for input A:
%      - logical, char, all numeric classes
%      - cell arrays of strings
%      -- 'rows' option is not supported for cell arrays
%      - objects with methods SORT (SORTROWS for the 'rows' option) and NE
%      -- including heterogeneous arrays
%
%   See also UNION, INTERSECT, SETDIFF, SETXOR, ISMEMBER, SORT, SORTROWS.

%   Copyright 1984-2012 The MathWorks, Inc.
%   $Revision: 1.24.4.13 $  $Date: 2012/07/28 23:09:26 $

% Determine the number of outputs requested.

if nargout == 0
    nlhs = 1;
else
    nlhs = nargout;
end

%narginchk(1,4);
nrhs = nargin;
if nrhs == 1
    [varargout{1:nlhs}] = uniqueR2012a(varargin{:});
else
    % acceptable combinations, with optional inputs denoted in []
    % unique(A, ['rows'], ['first'/'last'], ['legacy'/'R2012a']),
    % where the position of 'rows' and 'first'/'last' may be reversed
    % unique(A, ['rows'], ['sorted'/'stable']),
    % where the position of 'rows' and 'sorted'/'stable' may be reversed
    nflagvals = 7;
    flagvals = {'rows' 'first' 'last' 'sorted' 'stable' 'legacy' 'R2012a'};
    % When a flag is found, note the index into varargin where it was found
    flaginds = zeros(1,nflagvals);
    for i = 2:nrhs
        flag = varargin{i};
        foundflag = strcmpi(flag,flagvals);
        if ~any(foundflag)
            if ischar(flag)
                error(message('MATLAB:UNIQUE:UnknownFlag',flag));
            else
                error(message('MATLAB:UNIQUE:UnknownInput'));
            end
        end
        % Only 1 occurrence of each allowed flag value
        if flaginds(foundflag)
            error(message('MATLAB:UNIQUE:RepeatedFlag',flag));
        end
        flaginds(foundflag) = i;
    end
    
    % Only 1 of each of the paired flags
    if flaginds(2) && flaginds(3)
        error(message('MATLAB:UNIQUE:OccurrenceConflict'))
    end
    if flaginds(4) && flaginds(5)
        error(message('MATLAB:UNIQUE:SetOrderConflict'))
    end
    if flaginds(6) && flaginds(7)
        error(message('MATLAB:UNIQUE:BehaviorConflict'))
    end
    % 'legacy' and 'R2012a' flags must be trailing
    if flaginds(6) && flaginds(6)~=nrhs
        error(message('MATLAB:UNIQUE:LegacyTrailing'))
    end
    if flaginds(7) && flaginds(7)~=nrhs
        error(message('MATLAB:UNIQUE:R2012aTrailing'))
    end
    
    if flaginds(4) || flaginds(5) % 'stable'/'sorted' specified
        if flaginds(2) || flaginds(3) % does not combine with 'first'/'last'
            error(message('MATLAB:UNIQUE:SetOrderOccurrence'))
        end
        if flaginds(6) || flaginds(7) % does not combine with 'legacy'/'R2012a'
            error(message('MATLAB:UNIQUE:SetOrderBehavior'))
        end
        [varargout{1:nlhs}] = uniqueR2012a(varargin{1},logical(flaginds(1:5)));
    elseif flaginds(7) % trailing 'R2012a' specified
        [varargout{1:nlhs}] = uniqueR2012a(varargin{1},logical(flaginds(1:5)));
    elseif flaginds(6) % trailing 'legacy' specified
        [varargout{1:nlhs}] = uniquelegacy(varargin{1},logical(flaginds(1:3)));
    else % 'R2012a' (default behavior, to be changed to 'R2012a' in future)
        [varargout{1:nlhs}] = uniqueR2012a(varargin{1},logical(flaginds(1:5)));
    end
end
end

function [b,ndx,pos] = uniquelegacy(a,options)
% 'legacy' flag implementation 

%flagvals = {'rows' 'first' 'last'};
if nargin == 1
    byrow = false;
    order = 'last';
else
    byrow = (options(1) > 0);
    if options(2) > 0
        order = 'first';
    else % if options(3) > 0 || sum(options(2:3) == 0)
        order = 'last';
    end
end

rows = size(a,1);
cols = size(a,2);

rowvec = (rows == 1) && (cols > 1);

numelA = numel(a);
nOut = nargout;

if ~isa(a,'opaque')
    
    if ~byrow
        
        % Handle empty: no elements.
        
        if (numelA == 0)
            % Predefine b to be of the correct type.
            b = a([]);
            if max(size(a)) > 0
                b = reshape(b,0,1);
                ndx = zeros(0,1);
                pos = zeros(0,1);
            else
                ndx = [];
                pos = [];
            end
            return
            
        elseif (numelA == 1)
            % Scalar A: return the existing value of A.
            b = a; ndx = 1; pos = 1;
            return
            
            % General handling.
        else
            
            % Convert to columns
            a = a(:);
            
            % Convert to double array for purposes of faster sorting.
            % Additionally, UNIQUE calls DIFF, which requires double input.
            
            whichclass = class(a);
            isdouble = strcmp(whichclass,'double');
            
            if ~isdouble
                a = double(a);
            end
            
            % Sort if unsorted.  Only check this for long lists.
            
            checksortcut = 1000;
            
            if numelA <= checksortcut || ~(issorted(a))
                if nOut <= 1
                    b = sort(a);
                else
                    [b,ndx] = sort(a);
                end
            else
                b = a;
                if nOut > 1
                    ndx = (1:numelA)';  % If presorted, indices are 1,2,3,...
                end
            end
            
            % d indicates the location of non-matching entries.
            
            db = diff(b);
            
            % Since DIFF returns NaN in both Inf and NaN cases,
            % use slower method of detection if NaN's detected in DIFF(b).
            % After sort, Infs or NaNs will be at ends of list only.
            
            if (isnan(db(1)) || isnan(db(numelA-1)))
                d = b(1:numelA-1) ~= b(2:numelA);
            else
                d = db ~= 0;
            end
            
            if order(1) == 'l' % 'last'
                d(numelA,1) = true; % Final element is always a member of unique list.
            else % order == 'first'
                d = [true; d];      % First element is always a member of unique list.
            end
            
            b = b(d);         % Create unique list by indexing into sorted list.
            
            if nOut == 3
                if order(1) == 'l' % 'last'
                    pos = cumsum([1;full(d)]);  % Lists position, starting at 1.
                    pos(numelA+1) = [];         % Remove extra element introduced by d.
                else % order == 'first'
                    pos = cumsum(full(d));      % Lists position, starting at 1.
                end
                pos(ndx) = pos;             % Re-reference POS to indexing of SORT.
            end
            
            % Create indices if needed.
            if nOut > 1
                ndx = ndx(d);
            end
            
            % Re-convert to correct output data type using FEVAL.
            if ~isdouble
                b = feval(whichclass,b);
            end
        end
        
        % If row vector, return as row vector.
        if rowvec
            b = b.';
            if nOut > 1
                ndx = ndx.';
                if nOut > 2
                    pos = pos.';
                end
            end
        end
        
    else    % 'rows' case
        
        % Handle empty: no rows.
        
        if (rows == 0)
            % Predefine b to be of the correct type.
            b = a([]);
            ndx = [];
            pos = [];
            b = reshape(b,0,cols);
            if cols > 0
                ndx = reshape(ndx,0,1);
            end
            return
            
            % Handle scalar: one row.
            
        elseif (rows == 1)
            b = a; ndx = 1; pos = 1;
            return
        end
        
        % General handling.
        % Conversion to double not done: SORTROWS is slower for doubles
        % than other types.
        
        if nOut > 1
            [b,ndx] = sortrows(a);
        else
            b = sortrows(a);
        end
        
        % d indicates the location of non-matching entries.
        
        d = b(1:rows-1,:)~=b(2:rows,:);
        
        % d = 1 if differences between rows.  d = 0 if the rows are equal.
        
        d = any(d,2);
        if order(1) == 'l' % 'last'
            d(rows,1) = true; % Final row is always member of unique list.
        else % order = 'first'
            d = [true; d];    % First row is always a member of unique list.
        end
        
        b = b(d,:);         % Create unique list by indexing into sorted list.
        
        % Create position mapping vector using CUMSUM.
        
        if nOut == 3
            if order(1) == 'l' % 'last'
                pos = cumsum([1;full(d)]);  % Lists position, starting at 1.
                pos(rows+1) = [];           % Remove extra element introduced by d.
            else % order == 'first'
                pos = cumsum(full(d));      % Lists position, starting at 1.
            end
            pos(ndx) = pos;             % Re-reference POS to indexing of SORT.
        end
        
        % Create indices if needed.
        if nOut > 1
            ndx = ndx(d);
        end
    end
    
else
    % Handle objects that cannot be converted to doubles     
    if ~byrow
        
        % Handle empty: no elements.
        
        if (numelA == 0)
            % Predefine b to be of the correct type.
            b = a([]);
            if max(size(a)) > 0
                b = reshape(b,0,1);
                ndx = zeros(0,1);
                pos = zeros(0,1);
            else
                ndx = [];
                pos = [];
            end
            return
            
        elseif (numelA == 1)
            % Scalar A: return the existing value of A.
            b = a; ndx = 1; pos = 1;
            return
            
            % General handling.
        else
            
            % Convert to columns
            a = a(:);
            
            % Sort if unsorted.  Only check this for long lists.
            
            if nOut <= 1
                b = sort(a);
            else
                [b,ndx] = sort(a);
            end
            
            % d indicates the location of non-matching entries.
            
            d = b(1:numelA-1) ~= b(2:numelA);
            
            if order(1) == 'l' % 'last'
                d(numelA,1) = true; % Final element is always a member of unique list.
            else % order == 'first'
                d = [true; d];      % First element is always a member of unique list.
            end
            
            b = b(d);         % Create unique list by indexing into sorted list.
            
            if nOut == 3
                if order(1) == 'l' % 'last'
                    pos = cumsum([1;d]);        % Lists position, starting at 1.
                    pos(numelA+1) = [];         % Remove extra element introduced by d.
                else % order == 'first'
                    pos = cumsum(d);            % Lists position, starting at 1.
                end
                pos(ndx) = pos;             % Re-reference POS to indexing of SORT.
            end
            
            % Create indices if needed.
            if nOut > 1
                ndx = ndx(d);
            end
        end
        
        % If row vector, return as row vector.
        if rowvec
            b = b.';
            if nOut > 1
                ndx = ndx.';
                if nOut > 2
                    pos = pos.';
                end
            end
        end
        
    else    % 'rows' case
        
        % Handle empty: no rows.
        
        if (rows == 0)
            % Predefine b to be of the correct type.
            b = a([]);
            ndx = [];
            pos = [];
            b = reshape(b,0,cols);
            if cols > 0
                ndx = reshape(ndx,0,1);
            end
            return
            
            % Handle scalar: one row.
            
        elseif (rows == 1)
            b = a; ndx = 1; pos = 1;
            return
        end
        
        % General handling.
        
        if nOut > 1
            [b,ndx] = sortrows(a);
        else
            b = sortrows(a);
        end
        
        % d indicates the location of non-matching entries.
        
        d = b(1:rows-1,:)~=b(2:rows,:);
        
        % d = 1 if differences between rows.  d = 0 if the rows are equal.
        
        d = any(d,2);
        if order(1) == 'l' % 'last'
            d(rows,1) = true; % Final row is always member of unique list.
        else % order == 'first'
            d = [true; d];    % First row is always a member of unique list.
        end
        
        b = b(d,:);         % Create unique list by indexing into sorted list.
        
        % Create position mapping vector using CUMSUM.
        
        if nOut == 3
            pos = cumsum([1;d]);        % Lists position, starting at 1.
            pos(rows+1) = [];           % Remove extra element introduced by d.
            pos(ndx) = pos;             % Re-reference POS to indexing of SORT.
        end
        
        % Create indices if needed.
        if nOut > 1
            ndx = ndx(d);
        end
    end
end
end


function [c,indA,indC] = uniqueR2012a(a,options)
% 'R2012a' flag implementaion

% flagvals = {'rows' 'first' 'last' 'sorted' 'stable'};
if nargin == 1
    byrow = false;
    order = 'sorted';
else
    byrow = (options(1) > 0);
    if options(5) > 0
        order = 'stable';
    elseif options(3) > 0
        order = 'last';
    else % if options(4) > 0 || options(2) || sum(options(2:5) == 0)  
        order = 'sorted';   %'first' and 'sorted' do the same thing
    end
end

% Determine if A is a row vector.
rows = size(a,1);
cols = size(a,2);
rowvec = (rows == 1) && (cols > 1);

if ~byrow  % default case
    
    % Convert to column
    a = a(:);
    numelA = numel(a);
    
    % Sort A and get the indices if needed.
    if nargout > 1 || strcmp(order, 'stable')
        [sortA,indSortA] = sort(a);
    else
        sortA = sort(a);
    end
    
    % groupsSortA indicates the location of non-matching entries.
    if isnumeric(sortA) && (numelA > 1)
        dSortA = diff(sortA);
        if (isnan(dSortA(1)) || isnan(dSortA(numelA-1)))
            groupsSortA = sortA(1:numelA-1) ~= sortA(2:numelA);
        else
            groupsSortA = dSortA ~= 0;
        end
        
    else
        groupsSortA = sortA(1:numelA-1) ~= sortA(2:numelA);
    end
    
    if (numelA ~= 0)
        if strcmp(order, 'last') 
            groupsSortA = [groupsSortA; true];          % Final element is always a member of unique list.
        else  % if (strcmp(order, 'sorted') || strcmp(order, 'stable')) 
            groupsSortA = [true; groupsSortA];          % First element is always a member of unique list.
        end
    else
        groupsSortA = zeros(0,1);
    end
    
    % Extract unique elements.
    if strcmp(order, 'stable') 
        invIndSortA = indSortA;
        invIndSortA(invIndSortA) = 1:numelA;  % Find inverse permutation.
        logIndA = groupsSortA(invIndSortA);   % Create new logical by indexing into groupsSortA.
        c = a(logIndA);                       % Create unique list by indexing into unsorted a.
    else
        c = sortA(groupsSortA);         % Create unique list by indexing into sorted list.
    end
    
    % Find indA.
    if nargout > 1
        if strcmp(order, 'stable') 
            indA = find(logIndA);           % Find the indices of the unsorted logical.
        else
            indA = indSortA(groupsSortA);   % Find the indices of the sorted logical.
        end
    end
    
    % Find indC.
    if nargout == 3
        groupsSortA = full(groupsSortA);
        if numelA == 0
            indC = zeros(0,1);
        else
            switch order
                case 'last'
                    indC = cumsum([1;groupsSortA(1:end-1)]);            % Lists position, starting at 1.
                    indC(indSortA) = indC;                              % Re-reference indC to indexing of sortA.
                case 'sorted' 
                    indC = cumsum(groupsSortA);                         % Lists position, starting at 1.
                    indC(indSortA) = indC;                              % Re-reference indC to indexing of sortA.
                otherwise % 'stable'
                    [dummy,indSortC] = sort(c);                             % Sort C to get index.
                    
                    lengthGroupsSortA = diff(find([groupsSortA; true]));  % Determine how many of each of the above indices there are in IC.
                    
                    diffIndSortC = diff(indSortC);                        % Get the correct amount of each index.
                    diffIndSortC = [indSortC(1); diffIndSortC];
                    
                    indLengthGroupsSortA = cumsum([1; lengthGroupsSortA]);
                    indLengthGroupsSortA(end) = [];
                    
                    indCOrderedBySortA(indLengthGroupsSortA,1) = diffIndSortC;        % Since indCOrderedBySortA is not already established as a column,
                    if sum(lengthGroupsSortA) ~= length(indCOrderedBySortA);        % This is false if all the elements in A originally were unique and
                        indCOrderedBySortA(sum(lengthGroupsSortA),1) = 0;             % true if the original A had duplicates.
                    end
                    
                    indCOrderedBySortA = cumsum(indCOrderedBySortA);
                    indC = indCOrderedBySortA(invIndSortA);                 % Reorder the list of indices to the unsorted order.
            end
        end
    end
    
    % If A is row vector, return C as row vector.
    if rowvec
        c = c.';
    end
    
else    % 'rows' case
    if ~ismatrix(a)
        error(message('MATLAB:UNIQUE:ANotAMatrix'));
    end
    
    numRows = size(a,1);
    numCols = size(a,2);
    
    % Sort A and get the indices if needed.
    if nargout > 1 || strcmp(order, 'stable') 
        [sortA,indSortA] = sortrows(a);
    else
        sortA = sortrows(a);
    end
    
    % groupsSortA indicates the location of non-matching entries.
    groupsSortA = sortA(1:numRows-1,:) ~= sortA(2:numRows,:);
    groupsSortA = any(groupsSortA,2);
    if (numRows ~=0)
        if strcmp(order, 'last') 
            groupsSortA = [groupsSortA; true];          % Final row is always member of unique list.
        else  % if (strcmp(order, 'sorted') || strcmp(order, 'stable')) 
            groupsSortA = [true; groupsSortA];          % First row is always a member of unique list.
        end
    end
    
    % Extract Unique elements.
    if strcmp(order, 'stable')
        invIndSortA = indSortA;
        invIndSortA(invIndSortA) = 1:numRows;               % Find the inverse permutation of indSortA.
        logIndA = groupsSortA(invIndSortA);                 % Create new logical by indexing into groupsSortA.
        c = a(logIndA,:);                                   % Create unique list by indexing into unsorted a.
    else
        c = sortA(groupsSortA,:);         % Create unique list by indexing into sorted list.
    end
    
    % Find indA.
    if nargout > 1
        if strcmp(order, 'stable')
            indA = find(logIndA);           % Find the indices of the unsorted logical.
        else
            indA = indSortA(groupsSortA);   % Find the indices of the sorted logical.
        end
    end
    
    % Find indC.
    if nargout == 3
        groupsSortA = full(groupsSortA);
        switch order
            case 'last' 
                if (numRows == 0)
                    indC = cumsum(groupsSortA);                         % Empty A - use all of groupsSortA.
                    indC(indSortA) = indC;
                else
                    indC = cumsum([1;full(groupsSortA(1:end-1))]);      % Lists position, starting at 1.
                    indC(indSortA) = indC;                              % Re-reference indC to indexing of sortA.
                end
            case 'sorted' 
                indC = cumsum(groupsSortA);                             % Lists position, starting at 1.
                indC(indSortA) = indC;                                  % Re-reference indC to indexing of sortA.
            otherwise % 'stable'
                if numCols == 0
                    indC = ones(numRows,1);                             % For 'stable' ensure that empty A gives correct size and shape.
                elseif numRows == 0
                    indC = zeros(0,1);
                else
                    [dummy,indSortC] = sortrows(c);                         % Sort C to get index.
                    
                    lengthGroupsSortA = diff(find([groupsSortA; true]));    % Determine how many of each of the above indices there are in IC.
                    
                    diffIndSortC = diff(indSortC);
                    diffIndSortC = [indSortC(1); diffIndSortC];
                    
                    indLengthGroupsSortA = cumsum([1; lengthGroupsSortA]);  % Get the correct amount of each index.
                    indLengthGroupsSortA(end) = [];
                    
                    indCOrderedBySortA(indLengthGroupsSortA,1) = diffIndSortC;        % Since indCOrderedBySortA is not already established as a column,
                   
                    if sum(lengthGroupsSortA) ~= length(indCOrderedBySortA);
                        indCOrderedBySortA(sum(lengthGroupsSortA),1) = 0;
                    end
                    
                    indCOrderedBySortA = cumsum(indCOrderedBySortA);
                    indC = indCOrderedBySortA(invIndSortA);                 % Reorder the list of indices to the unsorted order.
                end
        end
    end
end
end
