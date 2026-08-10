function tests = test_temperature_processor_smoke
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.OriginalPath = path;
addpath(genpath(fullfile(repositoryRoot,'matlab')))
testCase.TestData.ExpectedFeatureNames = {
    'temperature_delta'
    'temperature_mean'
    'temperature_std'
    'temperature_delta_over_time'
    'temperature_linear_slope'
    'temperature_derivative_mean'
    'temperature_derivative_std'
    'temperature_derivative_linear_slope'
    };
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath)
end

function testConstructorDefaults(testCase)
processor = mpk.TemperatureProcessor();

verifyEqual(testCase,processor.ShowPlots,false)
verifyTrue(testCase,islogical(processor.ShowPlots))
end

function testProcessesOneSignalAndBuildsResult(testCase)
signal = localTemperatureSignal('temp',0,120);
processor = mpk.TemperatureProcessor();
result = processor.process(signal);

verifyClass(testCase,result,'mpk.ProcessingResult')
verifyEqual(testCase,result.Modality,'Temperature')
verifyEqual(testCase,result.ProcessorName,'mpk.TemperatureProcessor')
verifyEqual(testCase,result.SourceSignalLabel,signal.Label)
verifyEqual(testCase,result.Options,struct('ShowPlots',false))
verifyEqual(testCase,result.Metadata.LegacyFunction, ...
    'temperature_signal_processing')
verifyEqual(testCase,result.Metadata.SamplingFrequency, ...
    signal.SamplingFrequency)
verifyEqual(testCase,result.Metadata.NumberOfSamples,signal.count())
verifyEqual(testCase,result.Metadata.StartTimeSeconds,signal.Time(1))
verifyEqual(testCase,result.Metadata.EndTimeSeconds,signal.Time(end))
verifyTrue(testCase,iscell(result.Features))
verifySize(testCase,result.Features,[8 2])
verifyEqual(testCase,result.Features(:,1), ...
    testCase.TestData.ExpectedFeatureNames)
verifyTrue(testCase,all(cellfun( ...
    @(value) isnumeric(value) && isscalar(value),result.Features(:,2))))
verifyFalse(testCase,any(contains(result.Features(:,1),'__')))
end

function testRejectsInvalidSignalInputs(testCase)
processor = mpk.TemperatureProcessor();
recording = mpk.Recording('invalid_input');

verifyError(testCase,@() processor.process([1; 2; 3]), ...
    'mpk:TemperatureProcessor:InvalidSignal')
verifyError(testCase,@() processor.process(recording), ...
    'mpk:TemperatureProcessor:InvalidSignal')
end

function testRejectsInvalidShowPlots(testCase)
signal = localTemperatureSignal('temp',0,120);
processor = mpk.TemperatureProcessor();

processor.ShowPlots = 1;
verifyError(testCase,@() processor.process(signal), ...
    'mpk:TemperatureProcessor:InvalidShowPlots')

processor.ShowPlots = [true false];
verifyError(testCase,@() processor.process(signal), ...
    'mpk:TemperatureProcessor:InvalidShowPlots')
end

function testProceduralNumericalParity(testCase)
signal = localTemperatureSignal('temperature',0,120);
expected = temperature_signal_processing( ...
    signal.Samples, ...
    signal.Time, ...
    [1; signal.count()], ...
    signal.SamplingFrequency, ...
    false);

processor = mpk.TemperatureProcessor();
actual = processor.process(signal);

verifyEqual(testCase,actual.Features(:,1),expected(:,1))
verifyEqual(testCase,cell2mat(actual.Features(:,2)), ...
    cell2mat(expected(:,2)),'AbsTol',1e-12)
end

function testProcessingDoesNotMutateSignal(testCase)
metadata = struct('Sensor','synthetic','Channel',4);
signal = localTemperatureSignal('skin temperature',5,120,metadata);
originalSamples = signal.Samples;
originalTime = signal.Time;
originalSamplingFrequency = signal.SamplingFrequency;
originalLabel = signal.Label;
originalUnits = signal.Units;
originalMetadata = signal.Metadata;

processor = mpk.TemperatureProcessor();
processor.process(signal);

verifyEqual(testCase,signal.Samples,originalSamples)
verifyEqual(testCase,signal.Time,originalTime)
verifyEqual(testCase,signal.SamplingFrequency,originalSamplingFrequency)
verifyEqual(testCase,signal.Label,originalLabel)
verifyEqual(testCase,signal.Units,originalUnits)
verifyEqual(testCase,signal.Metadata,originalMetadata)
end

function testProcessorCanBeReusedForIndependentSignals(testCase)
processor = mpk.TemperatureProcessor();
firstSignal = localTemperatureSignal('first',0,120);
secondSignal = localTemperatureSignal('second',2,90);

firstResult = processor.process(firstSignal);
secondResult = processor.process(secondSignal);

verifyEqual(testCase,firstResult.SourceSignalLabel,'first')
verifyEqual(testCase,secondResult.SourceSignalLabel,'second')
verifyEqual(testCase,firstResult.Metadata.NumberOfSamples,1200)
verifyEqual(testCase,secondResult.Metadata.NumberOfSamples,900)
verifyNotEqual(testCase,firstResult.Features{2,2}, ...
    secondResult.Features{2,2})
end

function testRecordingProcessSignalIntegration(testCase)
signal = localTemperatureSignal('temperature channel',0,120);
recording = mpk.Recording('whole_signal');
recording.addSignal('temp',signal);
processor = mpk.TemperatureProcessor();

directResult = processor.process(signal);
recordingResult = recording.processSignal(processor,'temp');

localVerifyFeaturesEqual(testCase,recordingResult.Features, ...
    directResult.Features)
verifyEqual(testCase,recordingResult.SourceSignalLabel,signal.Label)
end

function testExtractedSignalCanBeProcessedDirectly(testCase)
recording = localPhasedRecording();
processor = mpk.TemperatureProcessor();
phaseSignal = recording.extractPhase('temp','segment_1');

result = processor.process(phaseSignal);

verifyClass(testCase,result,'mpk.ProcessingResult')
verifyEqual(testCase,result.SourceSignalLabel,phaseSignal.Label)
verifyEqual(testCase,result.Features(:,1), ...
    testCase.TestData.ExpectedFeatureNames)
end

function testRecordingProcessPhasesMatchesExplicitProcessing(testCase)
recording = localPhasedRecording();
processor = mpk.TemperatureProcessor();
originalSignal = recording.getSignal('temp');
originalSamples = originalSignal.Samples;
originalTime = originalSignal.Time;
originalIntervals = recording.PhaseIntervalsSeconds;
originalLabels = recording.PhaseLabels;

results = recording.processPhases(processor,'temp');

verifyClass(testCase,results,'containers.Map')
verifyEqual(testCase,sort(keys(results)),sort(originalLabels))
for index = 1:numel(originalLabels)
    phaseLabel = originalLabels{index};
    extracted = recording.extractPhase('temp',phaseLabel);
    directResult = processor.process(extracted);
    phaseResult = results(phaseLabel);
    verifyEqual(testCase,phaseResult.Features(:,1), ...
        testCase.TestData.ExpectedFeatureNames)
    verifyFalse(testCase,any(contains(phaseResult.Features(:,1),'__')))
    localVerifyFeaturesEqual(testCase,phaseResult.Features, ...
        directResult.Features)
end

verifyEqual(testCase,originalSignal.Samples,originalSamples)
verifyEqual(testCase,originalSignal.Time,originalTime)
verifyEqual(testCase,recording.PhaseIntervalsSeconds,originalIntervals)
verifyEqual(testCase,recording.PhaseLabels,originalLabels)
verifyNotEqual(testCase,results('segment_1').Features{2,2}, ...
    results('segment-2').Features{2,2})
end

function signal = localTemperatureSignal(label,startTime,durationSeconds,metadata)
if nargin < 4
    metadata = struct();
end
fs = 10;
nSamples = round(durationSeconds * fs);
relativeTime = (0:(nSamples-1))' / fs;
time = startTime + relativeTime;
samples = 36.2 + 0.006 * relativeTime ...
    + 0.04 * sin(2*pi*0.02*relativeTime);
signal = mpk.Signal(samples,time,fs,label,'degC',metadata);
end

function recording = localPhasedRecording()
signal = localTemperatureSignal('temperature channel',0,120);
recording = mpk.Recording('temperature_phases');
recording.addSignal('temp',signal);
recording.setPhases( ...
    [0, 60; 60, 120], ...
    {'segment_1','segment-2'});
end

function localVerifyFeaturesEqual(testCase,actual,expected)
verifyEqual(testCase,actual(:,1),expected(:,1))
verifyEqual(testCase,cell2mat(actual(:,2)),cell2mat(expected(:,2)), ...
    'AbsTol',1e-12)
end
