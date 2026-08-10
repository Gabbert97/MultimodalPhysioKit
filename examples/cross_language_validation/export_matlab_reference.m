function export_matlab_reference()
%EXPORT_MATLAB_REFERENCE Export final MATLAB features for notebook comparison.
%
% Run from the repository root after adding matlab/ and this example folder
% to the MATLAB path. EDA uses the separately installed external cvxEDA
% implementation available on the MATLAB path.

exampleDirectory = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(exampleDirectory));
addpath(genpath(fullfile(repositoryRoot,'matlab')))

fixturePath = fullfile(repositoryRoot,'examples','sample_data', ...
    'synthetic_multimodal_180s.hdf5');
outputDirectory = fullfile(repositoryRoot,'validation-output');
outputPath = fullfile(outputDirectory,'matlab_features.json');

recording = mpk.io.readBioSignalsPluxHdf5( ...
    fixturePath,'FnirsColumnOrder',[1 2]);
assert(recording.getSignal('fnirs_red').Metadata.ChannelId == 9, ...
    'CrossLanguageValidation:RedChannelMismatch', ...
    'The fNIRS red signal must come from channel 9.')
assert(recording.getSignal('fnirs_infrared').Metadata.ChannelId == 10, ...
    'CrossLanguageValidation:InfraredChannelMismatch', ...
    'The fNIRS infrared signal must come from channel 10.')

phaseIntervals = [0 60 120;60 120 180];
phaseLabels = {'phase_1','phase_2','phase_3'};
recording.setPhases(phaseIntervals,phaseLabels);
for index = 1:numel(phaseLabels)
    mapped = recording.phaseIndices('ecg',phaseLabels{index});
    assert(mapped(2)-mapped(1) == 6000, ...
        'CrossLanguageValidation:PhaseLengthMismatch', ...
        'Every phase must contain exactly 6000 samples.')
end

metadata = struct( ...
    'sampling_frequency',100, ...
    'phases',struct( ...
        'phase_1',[0 60], ...
        'phase_2',[60 120], ...
        'phase_3',[120 180]), ...
    'samples_per_phase',6000, ...
    'fnirs_age',25, ...
    'fnirs_red_channel',9, ...
    'fnirs_infrared_channel',10, ...
    'fnirs_baseline_phase','phase_1', ...
    'fnirs_baseline_interval_seconds',[10 50], ...
    'eda_decomposition','official_cvxEDA_v1.1.0', ...
    'minimum_scr_amplitude_us',0.01);

output = struct('metadata',metadata);

ecgProcessor = mpk.ECGProcessor('LabelFrequencyAnalysis',1);
output.ecg = localProcessSingleSignal( ...
    recording,ecgProcessor,'ecg',phaseLabels);

if exist('cvxEDA','file') ~= 2
    error('CrossLanguageValidation:MissingCvxEDA', ...
        ['The cross-language validation requires the external cvxEDA.m ' ...
         'implementation on the MATLAB path. Obtain cvxEDA separately, ' ...
         'review its GPL license, and add its directory to the MATLAB ' ...
         'path before running export_matlab_reference. See ' ...
         'docs/dependencies.md.'])
end
cvxedaFile = which('cvxEDA');
officialCvxedaSha256 = ...
    '66ad6b5dcd1d5131c2c3ad85e7dd15e4921a518027fbc84a097aa447f17fe179';
actualCvxedaSha256 = localSha256(cvxedaFile);
assert(strcmp(actualCvxedaSha256,officialCvxedaSha256), ...
    'CrossLanguageValidation:UnsupportedCvxEDA', ...
    ['The MATLAB validation requires the official, unmodified cvxEDA.m ' ...
     'from cvxEDA v1.1.0. The resolved file is %s.'],cvxedaFile)
output.metadata.cvxeda_matlab_sha256 = actualCvxedaSha256;
edaProcessor = mpk.EDAProcessor('MinimumScrAmplitude',0.01);
output.eda = localProcessSingleSignal( ...
    recording,edaProcessor,'eda',phaseLabels);
output.eda_accepted_event_counts = struct();
for index = 1:numel(phaseLabels)
    phaseName = phaseLabels{index};
    output.eda_accepted_event_counts.(phaseName) = round( ...
        output.eda.(phaseName).eda_scr_per_min * 59.99 / 60);
end

respirationProcessor = mpk.RespirationProcessor();
respirationProcessor.ShowPlots = false;
output.respiration = localProcessSingleSignal( ...
    recording,respirationProcessor,'rip',phaseLabels);

temperatureProcessor = mpk.TemperatureProcessor();
temperatureProcessor.ShowPlots = false;
output.temperature = localProcessSingleSignal( ...
    recording,temperatureProcessor,'temp',phaseLabels);

fnirsProcessor = mpk.FNIRSProcessor('AgeYears',25);
baselineRed = recording.extractPhase('fnirs_red','phase_1');
baselineInfrared = recording.extractPhase('fnirs_infrared','phase_1');
baseline = fnirsProcessor.computeBaselineValues( ...
    baselineRed,baselineInfrared);
assert(baseline.Diagnostics.EndIndex- ...
    baseline.Diagnostics.StartIndex+1 == 400, ...
    'CrossLanguageValidation:BaselineLengthMismatch', ...
    'The fNIRS baseline must contain exactly 400 samples at 10 Hz.')

fnirsResults = recording.processPhases( ...
    fnirsProcessor,{'fnirs_red','fnirs_infrared'}, ...
    'ProcessArguments',{ ...
        'BaselineRedValue',baseline.RedValue, ...
        'BaselineInfraredValue',baseline.InfraredValue});
output.fnirs = struct();
for index = 1:numel(phaseLabels)
    phaseName = phaseLabels{index};
    output.fnirs.(phaseName) = localFlattenFnirs( ...
        fnirsResults(phaseName).Features);
end

localAssertFinite(output)
if ~isfolder(outputDirectory)
    mkdir(outputDirectory)
end
fileIdentifier = fopen(outputPath,'w');
assert(fileIdentifier ~= -1, ...
    'CrossLanguageValidation:OutputOpenFailed', ...
    'Unable to open the MATLAB feature output file.')
cleanup = onCleanup(@() fclose(fileIdentifier));
fprintf(fileIdentifier,'%s\n',jsonencode(output,'PrettyPrint',true));
fprintf('MATLAB features written to %s\n',outputPath)
end

function exported = localProcessSingleSignal( ...
        recording,processor,signalName,phaseLabels)
results = recording.processPhases(processor,signalName);
exported = struct();
for index = 1:numel(phaseLabels)
    phaseName = phaseLabels{index};
    exported.(phaseName) = localFeatureTableToStruct( ...
        results(phaseName).Features);
end

end

function features = localFlattenFnirs(groups)
oxy = localFeatureTableToStruct(groups.Oxy);
deoxy = localFeatureTableToStruct(groups.Deoxy);
suffixes = fieldnames(oxy);
features = struct();
for index = 1:numel(suffixes)
    suffix = suffixes{index};
    features.(['fnirs_oxy_' suffix]) = oxy.(suffix);
    features.(['fnirs_deoxy_' suffix]) = deoxy.(suffix);
end
end

function features = localFeatureTableToStruct(table)
if iscell(table)
    names = string(table(:,1));
    values = zeros(numel(names),1);
    for index = 1:numel(names)
        value = table{index,2};
        if ischar(value) || isstring(value)
            value = str2double(value);
        end
        values(index) = double(value);
    end
elseif isstring(table)
    names = table(:,1);
    values = str2double(table(:,2));
else
    error('CrossLanguageValidation:UnsupportedFeatureTable', ...
        'Unsupported MATLAB feature table representation.')
end
features = struct();
for index = 1:numel(names)
    features.(char(names(index))) = values(index);
end
end

function localAssertFinite(value)
if isnumeric(value)
    assert(all(isfinite(double(value(:)))), ...
        'CrossLanguageValidation:NonfiniteOutput', ...
        'MATLAB feature output contains NaN or Inf.')
elseif isstruct(value)
    fields = fieldnames(value);
    for element = 1:numel(value)
        for index = 1:numel(fields)
            localAssertFinite(value(element).(fields{index}));
        end
    end
elseif iscell(value)
    for index = 1:numel(value)
        localAssertFinite(value{index});
    end
end
end

function digest = localSha256(filePath)
engine = javaMethod( ...
    'getInstance','java.security.MessageDigest','SHA-256');
fileIdentifier = fopen(filePath,'rb');
assert(fileIdentifier ~= -1, ...
    'CrossLanguageValidation:CvxedaReadFailed', ...
    'Unable to read the resolved cvxEDA.m file.')
cleanup = onCleanup(@() fclose(fileIdentifier));
bytes = fread(fileIdentifier,Inf,'*uint8');
engine.update(typecast(bytes(:),'int8'));
digestBytes = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(digestBytes,2).',1,[]));
end
