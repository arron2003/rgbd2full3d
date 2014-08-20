function dist = meanshift_dist2(X1, X2)

dist = sum((X1-X2).^2, 1);
