% Tests whether each element of one set is found in a second set.
%
% Args:
%   needles - the values we're searching for
%   haystack - the data that were searching in
%
% Test - a matrix the same size of 'needles' indicating whether each
%        element is found in 'haystack'.
function test = isin(needles, haystack)
  test = false(size(needles));
  for ii = 1 : length(haystack)
    test = test | (needles == haystack(ii));
  end
end
