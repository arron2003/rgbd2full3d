function h = sfigure(h)
% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%
% See also figure
error(nargchk(0, 1, nargin));

if nargin == 0
  h = figure;
  return
end

if ishandle(h)
  set(0, 'CurrentFigure', h);
else
  h = figure(h);
end