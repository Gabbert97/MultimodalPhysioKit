function tests = test_io_read_generic_hdf5
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

function testReadsExplicitlyMappedDatasets(testCase)
filePath = [tempname '.h5'];
cleanup = onCleanup(@() localDelete(filePath));

ecg = [0.1; 0.2; 0.3; 0.4];
ecgTime = [10; 10.5; 11; 11.5];
temperature = [36.5; 36.6; 36.7];
h5create(filePath,'/measurements/ecg',size(ecg));
h5write(filePath,'/measurements/ecg',ecg);
h5create(filePath,'/timing/ecg',size(ecgTime));
h5write(filePath,'/timing/ecg',ecgTime);
h5create(filePath,'/other/temperature',size(temperature));
h5write(filePath,'/other/temperature',temperature);

mapping(1) = struct('Name','ecg','DatasetPath','/measurements/ecg', ...
    'Units','mV','Label','lead_1','TimePath','/timing/ecg', ...
    'SamplingFrequency',[]);
mapping(2) = struct('Name','temperature', ...
    'DatasetPath','/other/temperature','Units','degC', ...
    'Label','skin_temperature','TimePath','', ...
    'SamplingFrequency',4);

recording = mpk.io.readGenericHdf5(filePath,mapping);

verifyClass(testCase,recording,'mpk.Recording')
verifyEqual(testCase,sort(recording.signalNames()), ...
    {'ecg'; 'temperature'})

ecgSignal = recording.getSignal('ecg');
verifyEqual(testCase,ecgSignal.Samples,ecg)
verifyEqual(testCase,ecgSignal.Time,ecgTime)
verifyEqual(testCase,ecgSignal.SamplingFrequency,2,'AbsTol',1e-12)
verifyEqual(testCase,ecgSignal.Units,'mV')
verifyEqual(testCase,ecgSignal.Label,'lead_1')

temperatureSignal = recording.getSignal('temperature');
verifyEqual(testCase,temperatureSignal.Samples,temperature)
verifyEqual(testCase,temperatureSignal.Time,[0; 0.25; 0.5])
verifyEqual(testCase,temperatureSignal.SamplingFrequency,4)
verifyEqual(testCase,temperatureSignal.Units,'degC')
verifyEqual(testCase,temperatureSignal.Label,'skin_temperature')
end

function testAcceptsMappingWithOnlySamplingFrequencyField(testCase)
filePath = [tempname '.h5'];
cleanup = onCleanup(@() localDelete(filePath));
samples = [1; 2; 3];
h5create(filePath,'/custom/signal',size(samples));
h5write(filePath,'/custom/signal',samples);
mapping = struct('Name','signal','DatasetPath','/custom/signal', ...
    'Units','a.u.','Label','custom','SamplingFrequency',5);

recording = mpk.io.readGenericHdf5(filePath,mapping);

signal = recording.getSignal('signal');
verifyEqual(testCase,signal.Samples,samples)
verifyEqual(testCase,signal.Time,[0; 0.2; 0.4],'AbsTol',1e-12)
end

function localDelete(filePath)
if isfile(filePath)
    delete(filePath)
end
end
