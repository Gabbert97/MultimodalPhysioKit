function data = eda_preprocess_signal(eda,time,fs,varargin)
%EDA_PREPROCESS_SIGNAL Validate, downsample, and optionally filter EDA.

parser = inputParser;
parser.FunctionName = 'eda_preprocess_signal';
addParameter(parser,'ApplyFilter',true, ...
    @(value) islogical(value) && isscalar(value));
parse(parser,varargin{:});

if ~isnumeric(eda) || ~isreal(eda) || ~isvector(eda) || ...
        ~isnumeric(time) || ~isreal(time) || ~isvector(time) || ...
        numel(eda) ~= numel(time)
    error('mpk:EDAProcessor:InvalidSignal', ...
        'EDA samples and timestamps must be equal-length numeric vectors.')
end
if ~isnumeric(fs) || ~isreal(fs) || ~isscalar(fs) || ~isfinite(fs)
    error('mpk:EDAProcessor:InvalidSamplingFrequency', ...
        'EDA SamplingFrequency must be an integer multiple of 100 Hz.')
end
ratio = fs / 100;
if ratio < 1 || abs(ratio - round(ratio)) > ...
        eps(max(1,abs(ratio))) * 8
    error('mpk:EDAProcessor:InvalidSamplingFrequency', ...
        'EDA SamplingFrequency must be an integer multiple of 100 Hz.')
end
if numel(eda) < 2
    error('mpk:EDAProcessor:SignalTooShort', ...
        'EDA signal is too short.')
end
if any(~isfinite(eda(:))) || any(~isfinite(time(:)))
    error('mpk:EDAProcessor:NonfiniteSignal', ...
        'EDA samples and timestamps must contain finite values.')
end
duration = time(end) - time(1);
if ~isfinite(duration) || duration <= 0
    error('mpk:EDAProcessor:NonpositiveDuration', ...
        'EDA signal duration must be positive.')
end

down = round(ratio);
indexes = 1:down:numel(eda);
edaDownsampled = antialias_downsample(double(eda(:)),down,100);
timeDownsampled = double(time(indexes));
edaFiltered = [];
if parser.Results.ApplyFilter
    try
        edaFiltered = filtering( ...
            edaDownsampled,100,0,0.5,0,5,0);
    catch exception
        errorWithCause = MException( ...
            'mpk:EDAProcessor:SignalTooShortForFiltering', ...
            'EDA signal is too short for low-pass filtering.');
        errorWithCause = addCause(errorWithCause,exception);
        throwAsCaller(errorWithCause)
    end
end

data = struct( ...
    'EdaDownsampled',edaDownsampled(:), ...
    'EdaFiltered',edaFiltered(:), ...
    'Time',timeDownsampled(:), ...
    'SamplingFrequency',100);
end
