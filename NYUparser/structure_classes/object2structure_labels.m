% Maps a given set of object labels to structure class labels.
%
% Args:
%   imgObjectLabels - HxW map of object labels where every pixel is in the
%                     range 1..C, where C is the number of classes or 0,
%                     indicating a missing label.
%   objectNames - Cx1 cell array of class names.
%
% Returns:
%   imgStructureLabels - HxW map of structure class labels where every
%                        pixel is in the range 0..4 where 1..4 indicate the
%                        structure class (see Consts.m) and 0 indicates a
%                        missing label.
function imgStructureLabels = object2structure_labels(imgObjectLabels, objectNames)
  obj2str = get_object_to_structure_class_map();
  
  % Grab all of the labels in the current label image.
  uniqueLabels = unique(imgObjectLabels);
  uniqueLabels(uniqueLabels == 0) = [];
  
  % Row, remap them all.
  imgStructureLabels = zeros(size(imgObjectLabels));
  for jj = 1 : numel(uniqueLabels)
    lblInd = uniqueLabels(jj);
    if ~isKey(obj2str, objectNames{lblInd})
      fprintf('\n');
      fprintf('Missing label for %s\n', objectNames{lblInd});
      fprintf('\n');
    else  
      imgStructureLabels(imgObjectLabels == lblInd) = obj2str(objectNames{lblInd});
    end
  end
end