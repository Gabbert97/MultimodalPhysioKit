function [features,intermediates] = eda_signal_processing( ...
        eda,time,fullSignalPhase,fs,varargin)
%EDA_SIGNAL_PROCESSING Process one complete EDA analysis window.
%
% [FEATURES,INTERMEDIATES] = EDA_SIGNAL_PROCESSING(EDA,TIME,[1;N],FS,...)
% downsamples and filters one complete window, resolves its cvxEDA
% decomposition, extracts eight scalar features, and optionally plots or
% saves decomposition diagnostics.

options = localParseOptions(varargin{:});
if ~isnumeric(fullSignalPhase) || ~isreal(fullSignalPhase) || ...
        ~isequal(size(fullSignalPhase),[2 1]) || ...
        ~isequal(fullSignalPhase,[1; numel(eda)])
    error('MultimodalPhysioKit:EDA:InvalidFullSignalPhase', ...
        'fullSignalPhase must equal [1; numel(eda)].')
end

preprocessed = eda_preprocess_signal(eda,time,fs);
decomposition = eda_resolve_decomposition( ...
    preprocessed,eda,time,fs,options);
features = eda_extract_features( ...
    decomposition,options.MinimumScrAmplitude);
eda_plot_decomposition(decomposition,options);

intermediates = struct( ...
    'eda_downsampled',preprocessed.EdaDownsampled, ...
    'eda_filtered',preprocessed.EdaFiltered, ...
    'time',preprocessed.Time, ...
    'sampling_frequency',preprocessed.SamplingFrequency, ...
    'scr',decomposition.Scr, ...
    'scl',decomposition.Scl, ...
    'residuals',decomposition.Residuals, ...
    'cvxeda_source',decomposition.CvxedaSource);
end

function options = localParseOptions(varargin)
parser = inputParser;
parser.FunctionName = 'eda_signal_processing';
addParameter(parser,'CvxedaPath','');
addParameter(parser,'DecompositionCallable',[]);
addParameter(parser,'CacheDirectory','');
addParameter(parser,'OutputLabel','eda');
addParameter(parser,'SignalLabel','eda');
addParameter(parser,'ShowPlots',false);
addParameter(parser,'SavePlots',false);
addParameter(parser,'OutputDirectory','');
addParameter(parser,'MinimumScrAmplitude',0.01);
parse(parser,varargin{:});

options = parser.Results;
options.CvxedaPath = localOptionalText( ...
    options.CvxedaPath,'CvxedaPath');
options.CacheDirectory = localOptionalText( ...
    options.CacheDirectory,'CacheDirectory');
options.OutputDirectory = localOptionalText( ...
    options.OutputDirectory,'OutputDirectory');
options.OutputLabel = localRequiredText( ...
    options.OutputLabel,'OutputLabel');
options.SignalLabel = localRequiredText( ...
    options.SignalLabel,'SignalLabel');
if ~(isempty(options.DecompositionCallable) || ...
        isa(options.DecompositionCallable,'function_handle'))
    error('mpk:EDAProcessor:InvalidDecompositionCallable', ...
        'DecompositionCallable must be empty or a function handle.')
end
if ~islogical(options.ShowPlots) || ~isscalar(options.ShowPlots)
    error('mpk:EDAProcessor:InvalidShowPlots', ...
        'ShowPlots must be a logical scalar.')
end
if ~islogical(options.SavePlots) || ~isscalar(options.SavePlots)
    error('mpk:EDAProcessor:InvalidSavePlots', ...
        'SavePlots must be a logical scalar.')
end
if options.SavePlots && isempty(options.OutputDirectory)
    error('mpk:EDAProcessor:MissingOutputDirectory', ...
        'OutputDirectory must be provided when SavePlots is true.')
end
if ~isnumeric(options.MinimumScrAmplitude) || ...
        ~isreal(options.MinimumScrAmplitude) || ...
        ~isscalar(options.MinimumScrAmplitude) || ...
        ~isfinite(options.MinimumScrAmplitude) || ...
        options.MinimumScrAmplitude < 0
    error('mpk:EDAProcessor:InvalidMinimumScrAmplitude', ...
        ['MinimumScrAmplitude must be a finite numeric scalar greater ' ...
         'than or equal to zero.'])
end
options.MinimumScrAmplitude = double(options.MinimumScrAmplitude);
end

function value = localOptionalText(rawValue,name)
if isempty(rawValue)
    value = '';
elseif ischar(rawValue) && isrow(rawValue)
    value = rawValue;
elseif isstring(rawValue) && isscalar(rawValue) && ~ismissing(rawValue)
    value = char(rawValue);
else
    error('mpk:EDAProcessor:InvalidConfiguration', ...
        '%s must be empty, a character vector, or a string scalar.',name)
end
end

function value = localRequiredText(rawValue,name)
if ischar(rawValue) && isrow(rawValue)
    value = rawValue;
elseif isstring(rawValue) && isscalar(rawValue) && ~ismissing(rawValue)
    value = char(rawValue);
else
    error('mpk:EDAProcessor:InvalidConfiguration', ...
        '%s must be a character vector or a string scalar.',name)
end
end
