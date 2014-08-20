function beta = Compute_beta(original_image,neighborSystem)

[M,N,S] = size(original_image);

count = 0;
beta = 0.0;
hm1 = original_image;
hm2 = (hm1(2:M,:,:)-hm1(1:(M-1),:,:)).^2;
[c1 c2 c3] = size(hm2);
count = count + (c1*c2*c3);
beta = beta + sum(sum(sum(hm2)));
hm2 = (hm1(:,2:N,:)-hm1(:,1:(N-1),:)).^2;
[c1 c2 c3] = size(hm2);
count = count + (c1*c2*c3);
beta = beta + sum(sum(sum(hm2)));
if (neighborSystem)
    hm2 = (hm1(2:M,2:N,:)-hm1(1:(M-1),1:(N-1),:)).^2;
    [c1 c2 c3] = size(hm2);
    count = count + (c1*c2*c3);
    beta = beta + sum(sum(sum(hm2)));
    hm2 = (hm1(1:(M-1),2:N,:)-hm1(2:M,1:(N-1),:)).^2;
    [c1 c2 c3] = size(hm2);
    count = count + (c1*c2*c3);
    beta = beta + sum(sum(sum(hm2)));
end

% final beta
beta_inv = 2 * (beta/count);
if (beta_inv == 0)
    beta_inv = epsilon;
end
beta = 1.0 / beta_inv;


end