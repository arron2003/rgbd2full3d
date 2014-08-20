function z = distSqr(x,y)

  if size(x,1)~=size(y,1), 
    error('size(x,1)~=size(y,1)'); 
  end

  [d,n] = size(x);
  [d,m] = size(y);

  z = x'*y;
  x2 = sum(x.^2,1)';
  y2 = sum(y.^2,1);
  for i = 1:m,
    z(:,i) = x2 + y2(i) - 2*z(:,i);
  end
end
