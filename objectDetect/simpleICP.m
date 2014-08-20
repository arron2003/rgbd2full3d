function M = simpleICP(p0, p1)
  % p0: static, p1: moving
  sqrErr = inf; convergence = 0.01; maxIter = 50;
  iter = 0; org_p1 = p1;
  while 1
    dist = pdist2(p0, p1); iter = iter+1;
    [mindist, minj] = min(dist, [], 2);
    newErr = sum(mindist);
    pm = p1(minj,:);
    M = mean(p0-pm); p1 = bsxfun(@plus, p1, M);
    if 0
      figure(2), clf; hold on; 
      scatter3(p0(:,1), p0(:,2), p0(:,3), 'r'); scatter3(p1(:,1), p1(:,2), p1(:,3), 'b');
      %line([p0(:,1) p1(:,1)], [p0(:,2) p1(:,2)], [p0(:,3) p1(:,3)]);
      fprintf('Iter %2d, %.3f \n', newErr);
    end
    
    if (iter>maxIter) || (newErr > sqrErr*(1-convergence))
      M=mean(p1-org_p1); break; 
    end
    sqrErr=newErr;
  end
  
end