function tests = test_io_recording_from_arrays
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.OriginalPath = path;
addpath(genpath(fullfile(repositoryRoot,'matlab')))
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath)
end

function testCreatesOneSignalFromSamplingFrequency(testCase)
samples = [1 2 3 4];
recording = mpk.io.recordingFromArrays(struct('ecg',samples), ...
    'SamplingFrequency',2);

verifyClass(testCase,recording,'mpk.Recording')
verifyTrue(testCase,recording.hasSignal('ecg'))
signal = recording.getSignal('ecg');
verifyClass(testCase,signal,'mpk.Signal')
verifyEqual(testCase,signal.Samples,samples(:))
verifyEqual(testCase,signal.Time,[0; 0.5; 1; 1.5])
verifyEqual(testCase,signal.SamplingFrequency,2)
verifyEqual(testCase,signal.Metadata, ...
    struct('SourceFormat','MATLAB arrays'))
end

function testCreatesMultipleSignalsFromExplicitTimes(testCase)
signals = struct('eda',[1 2 3],'temperature',[10 11 12 13]);
times = struct('eda',[5 5.5 6], ...
    'temperature',[0 0.25 0.5 0.75]);
recording = mpk.io.recordingFromArrays(signals,'Time',times);

verifyEqual(testCase,sort(recording.signalNames()), ...
    {'eda'; 'temperature'})
eda = recording.getSignal('eda');
temperature = recording.getSignal('temperature');
verifyEqual(testCase,eda.Time,times.eda(:))
verifyEqual(testCase,eda.SamplingFrequency,2,'AbsTol',1e-12)
verifyEqual(testCase,temperature.Time,times.temperature(:))
verifyEqual(testCase,temperature.SamplingFrequency,4,'AbsTol',1e-12)
verifyEqual(testCase,recording.Metadata, ...
    struct('SourceFormat','MATLAB arrays'))
end
