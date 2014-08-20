function [bestX, bestF] = beamSearch(unary, pairwise, complete, ios)
  n = numel(unary);
  nPool = 5;
  x0 = false(n, 1);
  % init
  fval = zeros(nPool*1000, 1);
  for i=1:(nPool*1000)
    [x{i} fval(i)] = hillClimb(unary, pairwise, complete, x0, 0, ~x0);
  end
  [~, maxi] = sort(fval);
  x = { x{maxi(1:nPool)} };
  fval = fval(maxi(1:nPool));
  
  % now try stuff
  maxIter = 100; cnt = 0;
  for iter=1:maxIter
    ridx = randperm(n);
    nochange = true;
    for obj=ridx
      for seed = 1:nPool
        x0 = x{seed};
        flag = false; newx = x0; newf = fval(seed);
        [tflag, tx, tf] = simpleSwap(unary, pairwise, complete, newx, newf, obj);
        if tflag, newx = tx; newf = tf; flag = true; end
        [tflag, tx, tf] = randSwap(unary, pairwise, complete, ios, newx, newf, obj);
        if tflag, newx = tx; newf = tf; flag = true; end

        if flag
          fprintf('%d = %f ~ %f\n', seed, fval(seed), newf);
          x{seed} = newx; fval(seed) = newf;
          nochange = false;
        end
      end
    end
    if nochange, break; end
  end
  
  [~, maxi] = sort(fval);
  x = { x{maxi(1:nPool)} };
  fval = fval(maxi(1:nPool));
  bestX = x{1}; bestF = fval(1);
  
  %[~, bestX, bestF] = randSwap(unary, pairwise, complete, ios, bestX, bestF, 88);
end

function [flag, x, f] = randSwap(unary, pairwise, complete, ios, x0, oldFval, act)
  x=x0; x(act)=true; 
  x( ios(act, :) )=false;
  newl = lossFunc(double(x), unary, pairwise, complete);
  [x, newl] = hillClimb(unary, pairwise, complete, x, newl, ~x);
  flag = newl < oldFval;
  f = newl;
end

function [flag, x, f] = simpleSwap(unary, pairwise, complete, x0, oldFval, act)
  x=x0; x(act)=true;
  newl = lossFunc(double(x), unary, pairwise, complete);
  flag = newl < oldFval;
  f = newl;
end

function [x, f] = hillClimb(unary, pairwise, complete, x0, oldFval, act)
  act = find(act);
  ridx = randperm(numel(act)); act = act(ridx)';
  f = oldFval;
  x = x0;
  tcomplete = max(bsxfun(@times, x0, complete));
  for i=act
    delta = unary(i) + 2*sum(pairwise(i, x)) - sum(max(complete(i,:)-tcomplete, 0));
    if delta<0, 
      x(i)=true; f = f + delta; 
      tcomplete = max(complete(i,:), tcomplete);
    end
  end
  c = tcomplete;
end

function l = lossFunc(x, unary, pairwise, complete)
  l = x'*unary'+x'*pairwise*x - sum( max(bsxfun(@times, x, complete)) );
end