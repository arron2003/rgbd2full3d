% Gets the time elapsed as a 3-vector.
%
% Args:
%   elapsedSeconds - the number of seconds elapsed
%
% Returns:
%   times - a 3-vector with the number of hours, minutes, and seconds.
function times = get_time_elapsed(elapsedSeconds)
  elapsedHours = floor(elapsedSeconds / 3600);
  
  elapsedSeconds = elapsedSeconds - elapsedHours * 3600;
  elapsedMinutes = floor(elapsedSeconds / 60);
  
  elapsedSeconds = elapsedSeconds - elapsedMinutes * 60;
  
  times = [elapsedHours elapsedMinutes elapsedSeconds];
end