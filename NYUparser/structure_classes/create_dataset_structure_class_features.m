% Creates a dataset of structure class features used for training the
% structure-class classifier.
%
% Args:
%   params - parameter structure.
function create_dataset_structure_class_features(params)
  Consts;

  % Load the train test split:
  load(consts.splitsPath, 'trainNdxs');

  %%
  trainData = [];
  trainLabels = [];

  testData = [];
  testLabels = [];

  for ii = 1 : consts.numImages
    fprintf('Loading structure features %d/%d\n', ii, consts.numImages);
    if ~consts.useImages(ii)
      continue;
    end

    params.seg.featureSet = 0;
    outFilename = sprintf(consts.structureFeaturesFilename, ...
      params.regionSrc, params.seg.featureSet, params.stage, ii);
    load(outFilename, 'regionFeatures');
    assert(~any(isnan(regionFeatures(:))));
    assert(~any(isinf(regionFeatures(:))));

    % Now, for each set of features, grab the appropriate labels.
    imgRegions = get_regions(ii, params);
    load(sprintf(consts.structureLabelsFilename, ii), 'imgStructureLabels');
    structureLabels = get_labels_from_regions(imgRegions, imgStructureLabels);
    assert(size(regionFeatures, 1) == numel(structureLabels));

    % Discard any data points where we're missing label information.
    inds = structureLabels > 0;
    regionFeatures = regionFeatures(inds, :);
    structureLabels = structureLabels(inds);

    if isin(ii, trainNdxs)
      trainData = [trainData; regionFeatures];
      trainLabels = [trainLabels; structureLabels];
    else
      testData = [testData; regionFeatures];
      testLabels = [testLabels; structureLabels];
    end
  end

  %%
  outFilename = sprintf(consts.structureFeaturesDataset, ...
      params.regionSrc, params.seg.featureSet, params.stage);
  fprintf('Saving file %s...', outFilename);
  save(outFilename, 'trainData', 'trainLabels', ...
      'testData', 'testLabels', '-v7.3');
  fprintf('DONE.\n');
end
