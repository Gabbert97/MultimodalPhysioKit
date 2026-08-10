function recording = readGenericHdf5(filePath,mapping)
%READGENERICHDF5 Read explicitly mapped signals from a generic HDF5 file.
%
% RECORDING = mpk.io.readGenericHdf5(FILEPATH,MAPPING) reads one signal for
% each element of the MAPPING struct array. No HDF5 group or attribute
% layout is assumed. Every mapping element must contain:
%   Name              - Name used to add the signal to the recording.
%   DatasetPath       - Absolute HDF5 path containing physical-unit data.
%   Units             - Signal units.
%   Label             - Descriptive signal label.
% and exactly one of these timing fields:
%   TimePath          - Absolute HDF5 path containing timestamps.
%   SamplingFrequency - Positive scalar sampling frequency.
%
% Exactly one of TimePath or SamplingFrequency must be supplied per signal.
% Dataset values are stored unchanged; no ADC or unit conversion is applied.

if ~(ischar(filePath) || (isstring(filePath) && isscalar(filePath)))
    error('mpk:io:readGenericHdf5:InvalidFilePath', ...
        'filePath must be a character vector or string scalar.')
end
filePath = char(filePath);
if ~isfile(filePath)
    error('mpk:io:readGenericHdf5:FileNotFound', ...
        'HDF5 file not found: %s',filePath)
end
if ~isstruct(mapping) || isempty(mapping)
    error('mpk:io:readGenericHdf5:InvalidMapping', ...
        'mapping must be a nonempty struct array.')
end

requiredFields = {'Name','DatasetPath','Units','Label'};
if ~all(isfield(mapping,requiredFields))
    error('mpk:io:readGenericHdf5:InvalidMapping', ...
        ['Each mapping entry must contain Name, DatasetPath, Units, ' ...
         'and Label.'])
end

recording = mpk.Recording('generic_hdf5_recording', ...
    struct('SourceFormat','Generic HDF5'));
for index = 1:numel(mapping)
    entry = mapping(index);
    name = localText(entry.Name,'Name',index,false);
    datasetPath = localText(entry.DatasetPath,'DatasetPath',index,false);
    units = localText(entry.Units,'Units',index,true);
    label = localText(entry.Label,'Label',index,true);

    hasTimePath = isfield(entry,'TimePath') && ...
        ~localIsMissing(entry.TimePath);
    hasSamplingFrequency = isfield(entry,'SamplingFrequency') && ...
        ~isempty(entry.SamplingFrequency);
    if hasTimePath == hasSamplingFrequency
        error('mpk:io:readGenericHdf5:InvalidTimingMapping', ...
            ['Mapping entry %d must supply exactly one of TimePath or ' ...
             'SamplingFrequency.'],index)
    end

    samples = h5read(filePath,datasetPath);
    if ~isnumeric(samples) || ~isvector(samples) || isempty(samples)
        error('mpk:io:readGenericHdf5:InvalidDataset', ...
            'DatasetPath for signal "%s" must contain a numeric vector.',name)
    end

    signalMetadata = struct('SourceFormat','Generic HDF5', ...
        'DatasetPath',datasetPath);
    if hasTimePath
        timePath = localText(entry.TimePath,'TimePath',index,false);
        time = h5read(filePath,timePath);
        [time,fs] = localValidateTime(time,samples,name);
        signalMetadata.TimePath = timePath;
    else
        fs = entry.SamplingFrequency;
        if ~isnumeric(fs) || ~isscalar(fs) || ~isfinite(fs) || fs <= 0
            error('mpk:io:readGenericHdf5:InvalidSamplingFrequency', ...
                ['SamplingFrequency for signal "%s" must be a positive ' ...
                 'finite numeric scalar.'],name)
        end
        time = (0:numel(samples)-1)' / fs;
    end

    signal = mpk.Signal(samples,time,fs,label,units,signalMetadata);
    recording.addSignal(name,signal);
end
end

function tf = localIsMissing(value)
tf = isempty(value) || ...
    (isstring(value) && isscalar(value) && strlength(value) == 0);
end

function value = localText(value,fieldName,index,allowEmpty)
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('mpk:io:readGenericHdf5:InvalidMapping', ...
        '%s in mapping entry %d must be text.',fieldName,index)
end
value = char(value);
if ~allowEmpty && isempty(value)
    error('mpk:io:readGenericHdf5:InvalidMapping', ...
        '%s in mapping entry %d must not be empty.',fieldName,index)
end
end

function [time,fs] = localValidateTime(time,samples,name)
if ~isnumeric(time) || ~isvector(time) || numel(time) ~= numel(samples)
    error('mpk:io:readGenericHdf5:InvalidTime', ...
        ['TimePath for signal "%s" must contain a numeric vector with ' ...
         'one timestamp per sample.'],name)
end
if numel(time) < 2 || any(~isfinite(time(:)))
    error('mpk:io:readGenericHdf5:InvalidTime', ...
        ['TimePath for signal "%s" must contain at least two finite ' ...
         'timestamps.'],name)
end

time = time(:);
intervals = diff(time);
if any(intervals <= 0)
    error('mpk:io:readGenericHdf5:InvalidTime', ...
        'Timestamps for signal "%s" must be strictly increasing.',name)
end
sampleInterval = median(intervals);
tolerance = max(10 * eps(max(abs(time))), ...
    1e-10 * abs(sampleInterval));
if any(abs(intervals - sampleInterval) > tolerance)
    error('mpk:io:readGenericHdf5:NonuniformTime', ...
        ['Timestamps for signal "%s" must be uniformly sampled so its ' ...
         'sampling frequency can be inferred.'],name)
end
fs = 1 / sampleInterval;
end
