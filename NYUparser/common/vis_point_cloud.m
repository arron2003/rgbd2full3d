function vis_point_cloud(points, colors)
  seq = randperm(size(points,1));
  numToPlot = min([5000, size(points,1)]);
  points2 = points(seq(1:numToPlot),:);

  switch size(points, 2)
    case 2
      
      X = points2(:,1);
      Y = points2(:,2);

      if nargin == 1
        scatter(X, Y, ones(numToPlot,1) * 30, 'r', 'filled');
      else
        if size(colors, 1) > 1
          colors = colors(seq(1:numToPlot), :);
        end
        scatter(X, Y, ones(numToPlot,1) * 30, colors, 'filled');
      end
      xlabel('x');
      ylabel('y');
      
    case 3
  
      X = points2(:,1);
      Y = points2(:,2);
      Z = points2(:,3);

      if nargin == 1
        scatter3(X, Y, Z, ones(numToPlot,1) * 30, 'r', 'filled');
      else
        if numel(colors) > 1
          colors = colors(seq(1:numToPlot), :);
        end
        scatter3(X, Y, Z, ones(numToPlot,1) * 30, 'filled', 'CData', colors);
        
%         scatter3(X, Y, Z, ones(numToPlot,1) * 30, 'filled', colors);
      end
      xlabel('x');
      ylabel('y');
      zlabel('z');
    otherwise 
      error('Points must be either 2 or 3d');
  end
  
  axis equal;
end
