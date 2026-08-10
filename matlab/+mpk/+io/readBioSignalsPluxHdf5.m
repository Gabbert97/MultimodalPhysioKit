function recording = readBioSignalsPluxHdf5(filePath,varargin)
%READBIOSIGNALSPLUXHDF5 Read sanitized BioSignalsPlux/OpenSignals HDF5 data.
%
% RECORDING = mpk.io.readBioSignalsPluxHdf5(FILEPATH) reads the first
% top-level recording group in FILEPATH and returns supported raw
% BioSignalsPlux channels as an mpk.Recording.
%
% RECORDING = mpk.io.readBioSignalsPluxHdf5(...,'FnirsColumnOrder',ORDER)
% selects which two hSpO2 channels become fnirs_red and fnirs_infrared.

if nargin < 1 || ~(ischar(filePath) || isstring(filePath)) || ...
        strlength(string(filePath)) == 0
    error('mpk:io:readBioSignalsPluxHdf5:InvalidFilePath', ...
        'filePath must be a nonempty character vector or string scalar.')
end
filePath = char(filePath);
if ~isfile(filePath)
    error('mpk:io:readBioSignalsPluxHdf5:FileNotFound', ...
        'HDF5 file not found.')
end

parser = inputParser;
parser.FunctionName = 'mpk.io.readBioSignalsPluxHdf5';
addParameter(parser,'FnirsColumnOrder',[1 2], ...
    @(value) isnumeric(value) && isvector(value) && numel(value) == 2 && ...
    all(isfinite(value)) && all(value == fix(value)) && all(value > 0));
parse(parser,varargin{:});
fnirsColumnOrder = reshape(double(parser.Results.FnirsColumnOrder),1,[]);

info = h5info(filePath);
if isempty(info.Groups)
    error('mpk:io:readBioSignalsPluxHdf5:InvalidHdf5', ...
        'The HDF5 file has no top-level recording group.')
end

datasetPath = info.Groups(1).Name;
channelIds = double(h5readatt(filePath,datasetPath,'channels'));
resolutionBits = double(h5readatt(filePath,datasetPath,'resolution'));
fs = localSamplingFrequency(filePath,datasetPath);

channelIds = channelIds(:);
resolutionBits = resolutionBits(:);
if isscalar(resolutionBits) && numel(channelIds) > 1
    resolutionBits = repmat(resolutionBits,numel(channelIds),1);
end
if numel(resolutionBits) ~= numel(channelIds)
    error('mpk:io:readBioSignalsPluxHdf5:InvalidChannelMetadata', ...
        'Channel and ADC resolution metadata have different lengths.')
end

recording = mpk.Recording('biosignalsplux',struct( ...
    'SourceFormat','BioSignalsPlux HDF5', ...
    'Dataset',datasetPath));

fnirsRaw = {};
fnirsChannelIds = [];
fnirsResolutionBits = [];
fnirsChannelPaths = {};
fnirsSensorLabels = {};

for index = 1:numel(channelIds)
    channelId = channelIds(index);
    channelPath = sprintf('%s/raw/channel_%d',datasetPath,channelId);
    sensorLabel = localSensorLabel(h5readatt(filePath,channelPath,'sensor'));
    canonicalSensor = localCanonicalSensor(sensorLabel);

    if strcmp(canonicalSensor,'hSpO2')
        fnirsRaw{end+1} = double(h5read(filePath,channelPath)); %#ok<AGROW>
        fnirsChannelIds(end+1,1) = channelId; %#ok<AGROW>
        fnirsResolutionBits(end+1,1) = resolutionBits(index); %#ok<AGROW>
        fnirsChannelPaths{end+1,1} = channelPath; %#ok<AGROW>
        fnirsSensorLabels{end+1,1} = sensorLabel; %#ok<AGROW>
        continue
    end

    signalName = localSignalName(canonicalSensor);
    if isempty(signalName) || recording.hasSignal(signalName)
        continue
    end

    raw = double(h5read(filePath,channelPath));
    raw = raw(:);
    converted = signal_conversion(raw,{canonicalSensor}, ...
        [channelId resolutionBits(index)]);
    samples = localConvertedSamples(converted,signalName);
    signal = mpk.Signal(samples,localTime(numel(samples),fs),fs, ...
        signalName,localUnits(signalName),localMetadata(datasetPath, ...
        channelId,channelPath,sensorLabel,resolutionBits(index)));
    recording.addSignal(signalName,signal);
end

if ~isempty(fnirsRaw)
    localAddFnirsSignals(recording,fnirsRaw,fnirsChannelIds, ...
        fnirsResolutionBits,fnirsChannelPaths,fnirsSensorLabels, ...
        fnirsColumnOrder,datasetPath,fs);
end

if isempty(recording.signalNames())
    error('mpk:io:readBioSignalsPluxHdf5:NoSupportedChannels', ...
        'No supported BioSignalsPlux channels were found in the recording group.')
end
end

function fs = localSamplingFrequency(filePath,datasetPath)
try
    fs = double(h5readatt(filePath,datasetPath,'sampling rate'));
catch
    fs = double(h5readatt(filePath,datasetPath,'fs'));
end
if ~isscalar(fs) || ~isfinite(fs) || fs <= 0
    error('mpk:io:readBioSignalsPluxHdf5:InvalidSamplingFrequency', ...
        'Sampling frequency must be a positive scalar.')
end
end

function label = localSensorLabel(value)
if iscell(value)
    value = value{1};
end
if isnumeric(value)
    label = char(value(:)');
else
    label = char(string(value));
end
label(label == char(0)) = [];
label = strtrim(label);
end

function canonicalSensor = localCanonicalSensor(sensorLabel)
if strcmpi(sensorLabel,'ECG')
    canonicalSensor = 'ECG';
elseif strcmpi(sensorLabel,'EDA')
    canonicalSensor = 'EDA';
elseif strcmpi(sensorLabel,'TEMP')
    canonicalSensor = 'TEMP';
elseif strcmpi(sensorLabel,'RIP')
    canonicalSensor = 'RIP';
elseif strcmpi(sensorLabel,'RESPIRATION')
    canonicalSensor = 'RESPIRATION';
elseif strcmpi(sensorLabel,'hSpO2')
    canonicalSensor = 'hSpO2';
else
    canonicalSensor = '';
end
end

function signalName = localSignalName(canonicalSensor)
switch canonicalSensor
    case 'ECG'
        signalName = 'ecg';
    case 'EDA'
        signalName = 'eda';
    case 'TEMP'
        signalName = 'temp';
    case 'RIP'
        signalName = 'rip';
    case 'RESPIRATION'
        signalName = 'respiration';
    otherwise
        signalName = '';
end
end

function units = localUnits(signalName)
switch signalName
    case 'ecg'
        units = 'mV';
    case 'eda'
        units = 'uS';
    case 'temp'
        units = 'degC';
    case 'rip'
        units = 'percent full scale';
    case 'respiration'
        units = 'V';
    case {'fnirs_red','fnirs_infrared'}
        units = 'uA';
    otherwise
        units = '';
end
end

function samples = localConvertedSamples(converted,signalName)
if ~isfield(converted,signalName)
    error('mpk:io:readBioSignalsPluxHdf5:ConversionFailed', ...
        'signal_conversion did not return the expected signal "%s".', ...
        signalName)
end
samples = converted.(signalName);
samples = samples(:);
end

function time = localTime(sampleCount,fs)
time = (0:(sampleCount - 1)).' ./ fs;
end

function metadata = localMetadata(datasetPath,channelId,channelPath, ...
        sensorLabel,resolutionBits)
metadata = struct( ...
    'SourceFormat','BioSignalsPlux HDF5', ...
    'Dataset',datasetPath, ...
    'ChannelId',channelId, ...
    'ChannelPath',channelPath, ...
    'SensorLabel',sensorLabel, ...
    'AdcResolutionBits',resolutionBits);
end

function localAddFnirsSignals(recording,fnirsRaw,channelIds,resolutionBits, ...
        channelPaths,sensorLabels,fnirsColumnOrder,datasetPath,fs)
if numel(fnirsRaw) < max(fnirsColumnOrder)
    error('mpk:io:readBioSignalsPluxHdf5:MissingFnirsChannels', ...
        'The requested FnirsColumnOrder requires more hSpO2 channels than the file contains.')
end

selectedRaw = fnirsRaw(fnirsColumnOrder);
sampleCounts = cellfun(@numel,selectedRaw);
if numel(unique(sampleCounts)) ~= 1
    error('mpk:io:readBioSignalsPluxHdf5:FnirsLengthMismatch', ...
        'Selected hSpO2 channels must have equal lengths.')
end

raw = zeros(sampleCounts(1),numel(selectedRaw));
for index = 1:numel(selectedRaw)
    raw(:,index) = selectedRaw{index}(:);
end

selectedChannelIds = channelIds(fnirsColumnOrder);
selectedResolutionBits = resolutionBits(fnirsColumnOrder);
converted = signal_conversion(raw,{'hSpO2','hSpO2'}, ...
    [selectedChannelIds selectedResolutionBits]);
if ~isfield(converted,'fnirs') || size(converted.fnirs,2) < 2
    error('mpk:io:readBioSignalsPluxHdf5:ConversionFailed', ...
        'signal_conversion did not return two fNIRS channels.')
end

names = {'fnirs_red','fnirs_infrared'};
for index = 1:2
    signalName = names{index};
    selectedIndex = fnirsColumnOrder(index);
    samples = converted.fnirs(:,index);
    metadata = localMetadata(datasetPath,selectedChannelIds(index), ...
        channelPaths{selectedIndex},sensorLabels{selectedIndex}, ...
        selectedResolutionBits(index));
    signal = mpk.Signal(samples,localTime(numel(samples),fs),fs, ...
        signalName,localUnits(signalName),metadata);
    recording.addSignal(signalName,signal);
end
end
