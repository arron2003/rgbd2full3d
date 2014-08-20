function parfor_save(outFilename, conf)
  f = fieldnames(conf);
  for i = 1:numel(f)
    eval([f{i} '=conf.' f{i} ';']);
  end
  save(outFilename, f{:});
end
