function l = lossFunc(x, unary, pairwise, complete)
  l = x'*unary'+x'*pairwise*x - sum( max(bsxfun(@times, x, complete)) );
end