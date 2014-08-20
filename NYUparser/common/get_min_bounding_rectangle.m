% Finds the minimum bounding rectangle using the method of rotating
% calipers.
%
% Args:
%   points2d - a 2D point cloud, an Nx2 matrix.
%   numTheta - the number of angles to check.
%
% Returns:
%   bb2d - struct containing the 2D basis and coefficients of the bounding
%          rectangle for the given 2D point cloud.
function bb2d = get_min_bounding_rectangle(points2d, numTheta)
  assert(ndims(points2d) == 2);
  assert(size(points2d,2) == 2);

  if nargin < 2
    numTheta = 5;
  end
  
  DEBUG = 0;

  N = size(points2d,1);

  % Pull the centroid before we shift everything to the origin.
  centroidOrig = mean(points2d);
  
  % Move the cloud to the mean.
  points2dCentered = points2d - repmat(centroidOrig, [N, 1]);
  
  % First, start by fitting the points using PCA.
  basis = princomp(points2dCentered);
  
  axis = basis(:,1);
  thetaAxis = atan2(axis(2), axis(1));
  
  % Show the prin axis.
  if DEBUG
    sfigure(1);
    vis_point_cloud(points2dCentered);
    vis_line(zeros(2,1), axis);
    title('Principle Axis');
    pause
  end
  
  % Now rotate the points so that the primary axis is aligned with the
  % positive X-axis.
  prinRot = [cos(-thetaAxis) -sin(-thetaAxis);
             sin(-thetaAxis)  cos(-thetaAxis);];
  points2dAligned = (prinRot * points2dCentered')';
  
  % Show the newly aligned point cloud.
  if DEBUG
    sfigure(2);
    vis_point_cloud(points2dAligned);
    title('Post alignment Axis');
    pause;
  end
  
  % Now that we have a starting point, find the bounding box with minimal
  % area.
  minArea = 10e10;
  rect = [];
  bestFitTheta = [];
  
  step = pi/(2*(numTheta-1));
  
  for theta = -pi/4 : step : pi/4
    
    % Rotate the points, not the basis. This will make it easier to fit a
    % rectangle to the points.
    R = [cos(theta) -sin(theta)
         sin(theta)  cos(theta)];
    points2dTmp = (R * points2dAligned')';
    
    % Now, find the axis aligned bounding box.
    minX = min(points2dTmp(:,1));
    maxX = max(points2dTmp(:,1));
    minY = min(points2dTmp(:,2));
    maxY = max(points2dTmp(:,2));
    
    if DEBUG
      sfigure(3);
      vis_point_cloud(points2dTmp);
      
      hold on;
      scatter(minX, minY, 30, 'b', 'filled');
      scatter(minX, maxY, 30, 'b', 'filled');
      scatter(maxX, minY, 30, 'b', 'filled');
      scatter(maxX, maxY, 30, 'b', 'filled');
      title(sprintf('Degrees=%f', 180*theta/pi));
      hold off;
      pause;
    end

    % Calculate the area.
    area = (maxY - minY) * (maxX - minX);
    if area < minArea
      if DEBUG
        sfigure(4);
        vis_point_cloud(points2dTmp);
        hold on;
        scatter(minX, minY, 30, 'b', 'filled');
        scatter(minX, maxY, 30, 'b', 'filled');
        scatter(maxX, minY, 30, 'b', 'filled');
        scatter(maxX, maxY, 30, 'b', 'filled');
        title('Best Fit');
        hold off;
      end

      rect = [minX maxX minY maxY];
      minArea = area;
      bestFitTheta = theta;
    end
  end
  
  %%
  bb2d = struct();
  bb2d.basis = [1 0; 0 1];
  bb2d.coeffs = [(rect(2) - rect(1))/2 (rect(4) - rect(3))/2];
  
  % Pretend the centroid is really at 0. This is necessary for the
  % rotation to be correct.
  bb2d.centroid = bb2d.coeffs + [rect(1), rect(3)];

  % Now, perform 2 rotations to the axes.
  theta2 = thetaAxis - bestFitTheta;
  
  R = [cos(theta2) -sin(theta2);
       sin(theta2)  cos(theta2)];
     
  yy = (R * bb2d.centroid')';

  bb2d.basis = (R * bb2d.basis')';
  bb2d.centroid = yy + centroidOrig;
end
