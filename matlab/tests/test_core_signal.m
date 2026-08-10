function tests = test_core_signal
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

function testConstructsSignalWithMetadata(testCase)
samples = [1 2 3];
time = [0 0.5 1];
metadata = struct('Source','unit_test');

signal = mpk.Signal(samples,time,2,'ecg','mV',metadata);

verifyClass(testCase,signal,'mpk.Signal')
verifyEqual(testCase,signal.Samples,[1; 2; 3])
verifyEqual(testCase,signal.Time,[0; 0.5; 1])
verifyEqual(testCase,signal.SamplingFrequency,2)
verifyEqual(testCase,signal.Label,'ecg')
verifyEqual(testCase,signal.Units,'mV')
verifyEqual(testCase,signal.Metadata,metadata)
verifyEqual(testCase,signal.count(),3)
verifyEqual(testCase,signal.durationSeconds(),1)
end

function testUsesDefaultOptionalFields(testCase)
signal = mpk.Signal([1; 2],[0; 1],1);

verifyEqual(testCase,signal.Label,'')
verifyEqual(testCase,signal.Units,'')
verifyEqual(testCase,signal.Metadata,struct())
end

function testRejectsInvalidSamples(testCase)
verifyError(testCase,@() mpk.Signal(ones(2),[0; 1],1), ...
    'mpk:Signal:InvalidSamples')
end

function testRejectsInvalidTime(testCase)
verifyError(testCase,@() mpk.Signal([1; 2],ones(2),1), ...
    'mpk:Signal:InvalidTime')
end

function testRejectsLengthMismatch(testCase)
verifyError(testCase,@() mpk.Signal([1; 2; 3],[0; 1],1), ...
    'mpk:Signal:LengthMismatch')
end

function testRejectsInvalidSamplingFrequency(testCase)
verifyError(testCase,@() mpk.Signal([1; 2],[0; 1],0), ...
    'mpk:Signal:InvalidSamplingFrequency')
verifyError(testCase,@() mpk.Signal([1; 2],[0; 1],Inf), ...
    'mpk:Signal:InvalidSamplingFrequency')
end

function testRejectsInvalidMetadata(testCase)
verifyError(testCase,@() mpk.Signal([1; 2],[0; 1],1,'ecg','mV','private'), ...
    'mpk:Signal:InvalidMetadata')
end
