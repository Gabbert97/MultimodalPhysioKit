function tests = test_io_read_biosignalsplux_hdf5
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.RepositoryRoot = repositoryRoot;
testCase.TestData.OriginalPath = path;
addpath(genpath(fullfile(repositoryRoot,'matlab')))
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath)
end

function testReadsTrackedSyntheticSample(testCase)
samplePath = fullfile(testCase.TestData.RepositoryRoot,'examples', ...
    'sample_data','synthetic_multimodal_180s.hdf5');
verifyTrue(testCase,isfile(samplePath), ...
    'The tracked synthetic sample fixture is missing.')

recording = mpk.io.readBioSignalsPluxHdf5(samplePath);

verifyClass(testCase,recording,'mpk.Recording')
verifyEqual(testCase,recording.Metadata.SourceFormat, ...
    'BioSignalsPlux HDF5')
verifyTrue(testCase,isfield(recording.Metadata,'Dataset'))
verifyFalse(testCase,isfield(recording.Metadata,'SourcePath'))
verifyFalse(testCase,isfield(recording.Metadata,'FilePath'))

verifyTrue(testCase,recording.hasSignal('ecg'))
verifyTrue(testCase,recording.hasSignal('eda'))
verifyTrue(testCase,recording.hasSignal('temp'))
verifyTrue(testCase,recording.hasSignal('fnirs_red'))
verifyTrue(testCase,recording.hasSignal('fnirs_infrared'))
verifyTrue(testCase,recording.hasSignal('rip') || ...
    recording.hasSignal('respiration'))

localVerifySignal(testCase,recording.getSignal('ecg'),'ecg','mV')
localVerifySignal(testCase,recording.getSignal('eda'),'eda','uS')
localVerifySignal(testCase,recording.getSignal('temp'),'temp','degC')
localVerifySignal(testCase,recording.getSignal('fnirs_red'), ...
    'fnirs_red','uA')
localVerifySignal(testCase,recording.getSignal('fnirs_infrared'), ...
    'fnirs_infrared','uA')
if recording.hasSignal('rip')
    localVerifySignal(testCase,recording.getSignal('rip'), ...
        'rip','percent full scale')
else
    localVerifySignal(testCase,recording.getSignal('respiration'), ...
        'respiration','V')
end

ecgSignal = recording.getSignal('ecg');
rawEcg = double(h5read(samplePath,ecgSignal.Metadata.ChannelPath));
expectedEcg = signal_conversion(rawEcg(:),{'ECG'}, ...
    [ecgSignal.Metadata.ChannelId ecgSignal.Metadata.AdcResolutionBits]);
verifyEqual(testCase,ecgSignal.Samples,expectedEcg.ecg(:), ...
    'AbsTol',1e-12)
end

function testUsesFnirsColumnOrder(testCase)
samplePath = fullfile(testCase.TestData.RepositoryRoot,'examples', ...
    'sample_data','synthetic_multimodal_180s.hdf5');
verifyTrue(testCase,isfile(samplePath), ...
    'The tracked synthetic sample fixture is missing.')

defaultRecording = mpk.io.readBioSignalsPluxHdf5(samplePath);
swappedRecording = mpk.io.readBioSignalsPluxHdf5(samplePath, ...
    'FnirsColumnOrder',[2 1]);

defaultRed = defaultRecording.getSignal('fnirs_red');
defaultInfrared = defaultRecording.getSignal('fnirs_infrared');
swappedRed = swappedRecording.getSignal('fnirs_red');
swappedInfrared = swappedRecording.getSignal('fnirs_infrared');

verifyEqual(testCase,swappedRed.Metadata.ChannelId, ...
    defaultInfrared.Metadata.ChannelId)
verifyEqual(testCase,swappedInfrared.Metadata.ChannelId, ...
    defaultRed.Metadata.ChannelId)
verifyEqual(testCase,swappedRed.Samples,defaultInfrared.Samples)
verifyEqual(testCase,swappedInfrared.Samples,defaultRed.Samples)
end

function localVerifySignal(testCase,signal,expectedLabel,expectedUnits)
verifyClass(testCase,signal,'mpk.Signal')
verifyEqual(testCase,signal.Label,expectedLabel)
verifyEqual(testCase,signal.Units,expectedUnits)
verifyGreaterThan(testCase,signal.count(),0)
verifyEqual(testCase,numel(signal.Samples),numel(signal.Time))
verifyGreaterThan(testCase,signal.SamplingFrequency,0)

requiredFields = {'SourceFormat','Dataset','ChannelId','ChannelPath', ...
    'SensorLabel','AdcResolutionBits'};
for index = 1:numel(requiredFields)
    verifyTrue(testCase,isfield(signal.Metadata,requiredFields{index}), ...
        sprintf('Missing metadata field %s.',requiredFields{index}))
end
verifyEqual(testCase,signal.Metadata.SourceFormat,'BioSignalsPlux HDF5')
verifyFalse(testCase,isfield(signal.Metadata,'SourcePath'))
verifyFalse(testCase,isfield(signal.Metadata,'FilePath'))
end
