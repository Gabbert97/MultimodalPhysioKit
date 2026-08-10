function decomposition = eda_resolve_decomposition( ...
        preprocessed,originalEda,originalTime,originalFs,options)
%EDA_RESOLVE_DECOMPOSITION Resolve, validate, and optionally cache cvxEDA.

expectedLength = numel(preprocessed.EdaFiltered);
cachePath = localCachePath( ...
    originalEda,originalTime,originalFs,options);
if ~isempty(cachePath) && isfile(cachePath)
    cached = localLoadCache(cachePath,expectedLength);
    scr = cached.scr;
    scl = cached.scl;
    residuals = cached.residuals;
else
    [scr,scl,residuals,source] = localRunCvxEda( ...
        preprocessed.EdaFiltered,options);
    scr = localValidateArray(scr,expectedLength,'Scr');
    scl = localValidateArray(scl,expectedLength,'Scl');
    residuals = localValidateArray( ...
        residuals,expectedLength,'Residuals');
    if ~isempty(cachePath)
        cacheDirectory = fileparts(cachePath);
        if ~isfolder(cacheDirectory)
            mkdir(cacheDirectory)
        end
        save(cachePath,'scr','scl','residuals')
    end
end
if ~exist('source','var')
    source = localDependencySource(options);
end

decomposition = struct( ...
    'Scr',scr(:), ...
    'Scl',scl(:), ...
    'Residuals',residuals(:), ...
    'EdaFiltered',preprocessed.EdaFiltered(:), ...
    'EdaDownsampled',preprocessed.EdaDownsampled(:), ...
    'Time',preprocessed.Time(:), ...
    'SamplingFrequency',preprocessed.SamplingFrequency, ...
    'CvxedaSource',source);
end

function [scr,scl,residuals,source] = localRunCvxEda(edaFiltered,options)
delta = 1 / 100;
if ~isempty(options.DecompositionCallable)
    callable = options.DecompositionCallable;
    source = 'callable';
elseif ~isempty(options.CvxedaPath)
    target = options.CvxedaPath;
    if isfolder(target)
        cvxedaFile = fullfile(target,'cvxEDA.m');
        cvxedaDirectory = target;
    elseif isfile(target)
        [cvxedaDirectory,name,extension] = fileparts(target);
        if ~strcmpi([name extension],'cvxEDA.m')
            error('mpk:EDAProcessor:InvalidCvxedaPath', ...
                'CvxedaPath file must identify cvxEDA.m.')
        end
        cvxedaFile = target;
    else
        error('mpk:EDAProcessor:InvalidCvxedaPath', ...
            'CvxedaPath must identify cvxEDA.m or its directory.')
    end
    if ~isfile(cvxedaFile)
        error('mpk:EDAProcessor:InvalidCvxedaPath', ...
            'cvxEDA.m was not found at the configured path.')
    end
    originalPath = path;
    cleanup = onCleanup(@() localRestorePath(originalPath));
    addpath(cvxedaDirectory,'-begin')
    rehash path
    callable = @cvxEDA;
    source = 'path';
else
    if exist('cvxEDA','file') ~= 2
        error('mpk:EDAProcessor:MissingCvxEDA', ...
            ['EDA decomposition requires external cvxEDA. Add cvxEDA.m ' ...
             'to the MATLAB path, provide CvxedaPath, or provide ' ...
             'DecompositionCallable.'])
    end
    callable = @cvxEDA;
    source = 'import';
end
if ~isempty(options.DecompositionCallable)
    [scr,~,scl,~,~,residuals,~] = callable(edaFiltered,delta);
else
    [scr,~,scl,~,~,residuals,~] = callable( ...
        edaFiltered,delta, ...
        2.0,0.7,10.0,8e-4,1e-2);
end
end

function source = localDependencySource(options)
if ~isempty(options.DecompositionCallable)
    source = 'callable';
elseif ~isempty(options.CvxedaPath)
    source = 'path';
else
    source = 'import';
end
end

function value = localValidateArray(rawValue,expectedLength,name)
if ~isnumeric(rawValue) || ~isvector(rawValue) || ...
        numel(rawValue) ~= expectedLength
    error('mpk:EDAProcessor:InvalidDecompositionOutput', ...
        '%s must have the same length as filtered EDA.',name)
end
value = double(rawValue(:));
if any(~isfinite(value))
    error('mpk:EDAProcessor:InvalidDecompositionOutput', ...
        '%s must contain finite numeric values.',name)
end
end

function cachePath = localCachePath(eda,time,fs,options)
if isempty(options.CacheDirectory)
    cachePath = '';
    return
end
digest = localSignalDigest(eda,time,fs);
cachePath = fullfile(options.CacheDirectory, ...
    sprintf('%s_%s_decomposition.mat', ...
        options.OutputLabel,digest(1:16)));
end

function cached = localLoadCache(cachePath,expectedLength)
try
    loaded = load(cachePath);
catch exception
    errorWithCause = MException( ...
        'mpk:EDAProcessor:InvalidCache', ...
        'EDA decomposition cache could not be loaded.');
    errorWithCause = addCause(errorWithCause,exception);
    throwAsCaller(errorWithCause)
end
required = {'scr','scl','residuals'};
if ~all(isfield(loaded,required))
    error('mpk:EDAProcessor:InvalidCache', ...
        'EDA decomposition cache is missing required arrays.')
end
cached = struct( ...
    'scr',localValidateCacheArray( ...
        loaded.scr,expectedLength,'scr'), ...
    'scl',localValidateCacheArray( ...
        loaded.scl,expectedLength,'scl'), ...
    'residuals',localValidateCacheArray( ...
        loaded.residuals,expectedLength,'residuals'));
end

function value = localValidateCacheArray(rawValue,expectedLength,name)
if ~isnumeric(rawValue) || ~isvector(rawValue) || ...
        numel(rawValue) ~= expectedLength || ...
        any(~isfinite(rawValue(:)))
    error('mpk:EDAProcessor:InvalidCache', ...
        ['Cached %s must be a finite numeric vector with the ' ...
         'filtered EDA length.'],name)
end
value = double(rawValue(:));
end

function digest = localSignalDigest(eda,time,fs)
engine = javaMethod( ...
    'getInstance','java.security.MessageDigest','SHA-256');
localUpdateDigest(engine,double(eda(:)));
localUpdateDigest(engine,double(time(:)));
localUpdateDigest(engine,double(fs));
digestBytes = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(digestBytes,2).',1,[]));
end

function localUpdateDigest(engine,values)
bytes = typecast(values(:),'int8');
engine.update(bytes);
end

function localRestorePath(originalPath)
path(originalPath)
rehash path
end
