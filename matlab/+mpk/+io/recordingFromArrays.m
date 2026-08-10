function recording = recordingFromArrays(namedSignals,varargin)
%RECORDINGFROMARRAYS Create an mpk.Recording from named numeric arrays.
%
% RECORDING = mpk.io.recordingFromArrays(NAMEDSIGNALS,'SamplingFrequency',FS)
% creates time vectors starting at zero. NAMEDSIGNALS is a scalar struct
% whose fields contain numeric signal vectors. FS may be a positive scalar
% shared by every signal or a struct with one value per signal.
%
% RECORDING = mpk.io.recordingFromArrays(NAMEDSIGNALS,'Time',TIME) uses
% explicit, uniformly sampled time vectors and infers each sampling
% frequency. TIME may be a numeric vector shared by every signal or a
% struct with one vector per signal.
%
% Exactly one of Time or SamplingFrequency must be supplied. The resulting
% recording and signals contain only array-source metadata; no file paths or
% file-format state are retained.

if ~isstruct(namedSignals) || ~isscalar(namedSignals)
    error('mpk:io:recordingFromArrays:InvalidSignals', ...
        'namedSignals must be a scalar struct of numeric vectors.')
end
signalNames = fieldnames(namedSignals);
if isempty(signalNames)
    error('mpk:io:recordingFromArrays:NoSignals', ...
        'namedSignals must contain at least one signal.')
end

parser = inputParser;
parser.FunctionName = 'mpk.io.recordingFromArrays';
addParameter(parser,'Time',[], ...
    @(value) isempty(value) || isnumeric(value) || isstruct(value));
addParameter(parser,'SamplingFrequency',[], ...
    @(value) isempty(value) || isnumeric(value) || isstruct(value));
parse(parser,varargin{:});
options = parser.Results;

hasTime = ~isempty(options.Time);
hasSamplingFrequency = ~isempty(options.SamplingFrequency);
if hasTime == hasSamplingFrequency
    error('mpk:io:recordingFromArrays:TimingRequired', ...
        'Supply exactly one of Time or SamplingFrequency.')
end

sourceMetadata = struct('SourceFormat','MATLAB arrays');
recording = mpk.Recording('array_recording',sourceMetadata);

for index = 1:numel(signalNames)
    name = signalNames{index};
    samples = namedSignals.(name);
    if ~isnumeric(samples) || ~isvector(samples) || isempty(samples)
        error('mpk:io:recordingFromArrays:InvalidSignal', ...
            'Signal "%s" must be a nonempty numeric vector.',name)
    end

    if hasTime
        time = localValueForSignal(options.Time,name,'Time');
        [time,fs] = localValidateTime(time,samples,name);
    else
        fs = localValueForSignal(options.SamplingFrequency,name, ...
            'SamplingFrequency');
        if ~isnumeric(fs) || ~isscalar(fs) || ~isfinite(fs) || fs <= 0
            error('mpk:io:recordingFromArrays:InvalidSamplingFrequency', ...
                ['SamplingFrequency for signal "%s" must be a positive ' ...
                 'finite numeric scalar.'],name)
        end
        time = (0:numel(samples)-1)' / fs;
    end

    signal = mpk.Signal(samples,time,fs,name,'',sourceMetadata);
    recording.addSignal(name,signal);
end
end

function value = localValueForSignal(optionValue,name,optionName)
if isstruct(optionValue)
    if ~isscalar(optionValue) || ~isfield(optionValue,name)
        error('mpk:io:recordingFromArrays:MissingTiming', ...
            '%s does not contain a value for signal "%s".',optionName,name)
    end
    value = optionValue.(name);
else
    value = optionValue;
end
end

function [time,fs] = localValidateTime(time,samples,name)
if ~isnumeric(time) || ~isvector(time) || numel(time) ~= numel(samples)
    error('mpk:io:recordingFromArrays:InvalidTime', ...
        ['Time for signal "%s" must be a numeric vector with one value ' ...
         'per sample.'],name)
end
if numel(time) < 2 || any(~isfinite(time(:)))
    error('mpk:io:recordingFromArrays:InvalidTime', ...
        ['Time for signal "%s" must contain at least two finite ' ...
         'values.'],name)
end

time = time(:);
intervals = diff(time);
if any(intervals <= 0)
    error('mpk:io:recordingFromArrays:InvalidTime', ...
        'Time for signal "%s" must be strictly increasing.',name)
end
sampleInterval = median(intervals);
tolerance = max(10 * eps(max(abs(time))), ...
    1e-10 * abs(sampleInterval));
if any(abs(intervals - sampleInterval) > tolerance)
    error('mpk:io:recordingFromArrays:NonuniformTime', ...
        ['Time for signal "%s" must be uniformly sampled so its ' ...
         'sampling frequency can be inferred.'],name)
end
fs = 1 / sampleInterval;
end
