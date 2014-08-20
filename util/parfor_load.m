function varargout = parfor_load(inFilename, varargin)
  if nargout~=numel(varargin) || numel(varargin)==0
    error('Argument Number mismatch!');
  end
  t = load(inFilename, varargin{:});
  for i=1:numel(varargin)
    varargout{i} = getfield(t, varargin{i});
  end
end
