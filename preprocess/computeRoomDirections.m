function R = computeRoomDirections(data)
  % Credit: part of this is inspired by ITQ.m by Yunchao Gong
  R = eye(3);
  iter = 15;
  V = reshape(data.normal, [size(data.normal, 1)*size(data.normal, 2), 3]);
  idx = data.rawDepths~=0;
  V(~idx,:)=[];
  
  for ii = 1:iter
    nidx = sum(abs(V*R)>(1-2^(-ii/4)), 2)>=1;
    if sum(nidx)>1e4
      V=V(nidx,:);
    end
    sZ = V*R;
    Z = abs(sZ);
    UX = [Z(:,1)>Z(:,2)&Z(:,1)>Z(:,3), Z(:,2)>=Z(:,1)&Z(:,2)>=Z(:,3), Z(:,3)>=Z(:,1)&Z(:,3)>Z(:,2)];
    UX = UX .* sign(sZ);
    C = UX' * V;
    [UB,sigma,UA] = svd(C);
    R = UA * UB';
  end
end