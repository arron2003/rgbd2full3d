function w = nn_vectorize_weights(weights, biases)
  % Vectorize from the bottom and vectorize weights and then biases.
  w = weights{1}(:);
  for ii = 2 : length(weights)
    w = [w; weights{ii}(:)];
  end
  
  for ii = 1 : length(biases)
    w = [w; biases{ii}(:)];
  end
end
