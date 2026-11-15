function tests = test_respiration_processor_smoke
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.OriginalPath = path;
addpath(genpath(fullfile(repositoryRoot,'matlab')))
testCase.TestData.ExpectedFeatureNames = {
    'respiration_breath_rate_mean'
    'respiration_breath_rate_std'
    'respiration_inspiratory_time_mean'
    'respiration_inspiratory_time_std'
    'respiration_expiratory_time_mean'
    'respiration_expiratory_time_std'
    'respiration_timing_ratio_mean'
    'respiration_timing_ratio_std'
    'respiration_amplitude_mean'
    'respiration_amplitude_std'
    'respiration_minute_ventilation_mean'
    'respiration_minute_ventilation_std'
    'respiration_coefficient_of_variation'
    'respiration_autocorrelation'
    };
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath)
end

function testConstructorDefaults(testCase)
processor = mpk.RespirationProcessor();

verifyEqual(testCase,processor.ShowPlots,false)
verifyTrue(testCase,islogical(processor.ShowPlots))
end

function testProcessesOneSignalAndBuildsResult(testCase)
signal = localRespirationSignal('respiration',0,120,100);
processor = mpk.RespirationProcessor();
result = processor.process(signal);

verifyClass(testCase,result,'mpk.ProcessingResult')
verifyEqual(testCase,result.Modality,'Respiration')
verifyEqual(testCase,result.ProcessorName,'mpk.RespirationProcessor')
verifyEqual(testCase,result.SourceSignalLabel,signal.Label)
verifyEqual(testCase,result.Options,struct('ShowPlots',false))
verifyEqual(testCase,result.Metadata.LegacyFunction, ...
    'respiration_signal_processing')
verifyEqual(testCase,result.Metadata.OriginalSamplingFrequency,100)
verifyEqual(testCase,result.Metadata.ProcessedSamplingFrequency,100)
verifyEqual(testCase,result.Metadata.OriginalNumberOfSamples,12000)
verifyEqual(testCase,result.Metadata.ProcessedNumberOfSamples,12000)
verifyEqual(testCase,result.Metadata.StartTimeSeconds,signal.Time(1))
verifyEqual(testCase,result.Metadata.EndTimeSeconds,signal.Time(end))
verifyTrue(testCase,iscell(result.Features))
verifySize(testCase,result.Features,[14 2])
verifyEqual(testCase,result.Features(:,1), ...
    testCase.TestData.ExpectedFeatureNames)
verifyTrue(testCase,all(cellfun( ...
    @(value) isnumeric(value) && isscalar(value),result.Features(:,2))))
verifyFalse(testCase,any(contains(result.Features(:,1),'__')))
end

function testFeaturesUnaffectedByInterferenceAboveNewNyquist(testCase)
% 100.1 Hz interference would alias to 0.1 Hz under plain decimation from
% 1000 Hz to 100 Hz and distort breath detection and amplitudes.
fs = 1000;
time = (0:(120*fs-1))' / fs;
breathing = sin(2*pi*(0.25*time + 0.3*sin(2*pi*0.013*time)));
clean = mpk.Signal(breathing,time,fs,'rip','');
noisy = mpk.Signal(breathing + 0.5*sin(2*pi*100.1*time),time,fs,'rip','');
processor = mpk.RespirationProcessor();
cleanFeatures = processor.process(clean).Features;
noisyFeatures = processor.process(noisy).Features;
names = {'respiration_breath_rate_mean','respiration_amplitude_mean', ...
    'respiration_inspiratory_time_mean','respiration_expiratory_time_mean'};
for index = 1:numel(names)
    row = strcmp(cleanFeatures(:,1),names{index});
    verifyEqual(testCase,noisyFeatures{row,2},cleanFeatures{row,2}, ...
        'RelTol',1e-3,names{index})
end
end

function testProcessedSampleMetadataReflectsDownsampling(testCase)
signal = localRespirationSignal('belt',0,60,200);
processor = mpk.RespirationProcessor();
result = processor.process(signal);

verifyEqual(testCase,result.Metadata.OriginalSamplingFrequency,200)
verifyEqual(testCase,result.Metadata.ProcessedSamplingFrequency,100)
verifyEqual(testCase,result.Metadata.OriginalNumberOfSamples,12000)
verifyEqual(testCase,result.Metadata.ProcessedNumberOfSamples,6000)
end

function testRejectsNonSignalInputs(testCase)
processor = mpk.RespirationProcessor();
recording = mpk.Recording('invalid_input');

verifyError(testCase,@() processor.process([1; 2; 3]), ...
    'mpk:RespirationProcessor:InvalidSignal')
verifyError(testCase,@() processor.process(recording), ...
    'mpk:RespirationProcessor:InvalidSignal')
verifyError(testCase,@() processor.process(struct()), ...
    'mpk:RespirationProcessor:InvalidSignal')
end

function testRejectsInvalidShowPlots(testCase)
signal = localRespirationSignal('respiration',0,60,100);
processor = mpk.RespirationProcessor();

processor.ShowPlots = 1;
verifyError(testCase,@() processor.process(signal), ...
    'mpk:RespirationProcessor:InvalidShowPlots')

processor.ShowPlots = [true false];
verifyError(testCase,@() processor.process(signal), ...
    'mpk:RespirationProcessor:InvalidShowPlots')
end

function testProceduralSamplingFrequencyErrorPropagates(testCase)
signal = localRespirationSignal('respiration',0,60,150);
processor = mpk.RespirationProcessor();

verifyError(testCase,@() processor.process(signal), ...
    'MultimodalPhysioKit:InvalidSamplingFrequency')
end

function testProceduralNumericalParity(testCase)
signal = localRespirationSignal('respiration',0,120,100);
expected = respiration_signal_processing( ...
    signal.Samples, ...
    signal.Time, ...
    [1; signal.count()], ...
    signal.SamplingFrequency, ...
    false);

processor = mpk.RespirationProcessor();
actual = processor.process(signal);

localVerifyFeaturesEqual(testCase,actual.Features,expected)
end

function testProcessingDoesNotMutateSignal(testCase)
metadata = struct('Sensor','synthetic','Channel',2);
signal = localRespirationSignal('chest_belt',5,120,100,metadata);
originalSamples = signal.Samples;
originalTime = signal.Time;
originalSamplingFrequency = signal.SamplingFrequency;
originalLabel = signal.Label;
originalUnits = signal.Units;
originalMetadata = signal.Metadata;

processor = mpk.RespirationProcessor();
processor.process(signal);

verifyEqual(testCase,signal.Samples,originalSamples)
verifyEqual(testCase,signal.Time,originalTime)
verifyEqual(testCase,signal.SamplingFrequency,originalSamplingFrequency)
verifyEqual(testCase,signal.Label,originalLabel)
verifyEqual(testCase,signal.Units,originalUnits)
verifyEqual(testCase,signal.Metadata,originalMetadata)
end

function testProcessorIsIndependentOfSignalLabel(testCase)
processor = mpk.RespirationProcessor();
labels = {'respiration','rip','user_defined_belt'};

for index = 1:numel(labels)
    signal = localRespirationSignal(labels{index},0,60,100);
    result = processor.process(signal);
    verifyEqual(testCase,result.SourceSignalLabel,labels{index})
    verifyEqual(testCase,result.Features(:,1), ...
        testCase.TestData.ExpectedFeatureNames)
end
end

function testProcessorCanBeReusedForIndependentSignals(testCase)
processor = mpk.RespirationProcessor();
firstSignal = localRespirationSignal('first',0,60,100);
secondSignal = localRespirationSignal('second',3,90,100);

firstResult = processor.process(firstSignal);
secondResult = processor.process(secondSignal);

verifyEqual(testCase,firstResult.SourceSignalLabel,'first')
verifyEqual(testCase,secondResult.SourceSignalLabel,'second')
verifyEqual(testCase,firstResult.Metadata.OriginalNumberOfSamples,6000)
verifyEqual(testCase,secondResult.Metadata.OriginalNumberOfSamples,9000)
verifySize(testCase,firstResult.Features,[14 2])
verifySize(testCase,secondResult.Features,[14 2])
end

function testRecordingProcessSignalSupportsAnyStorageName(testCase)
signal = localRespirationSignal('chest_belt',0,120,100);
processor = mpk.RespirationProcessor();

respirationRecording = mpk.Recording('respiration_name');
respirationRecording.addSignal('respiration',signal);
respirationResult = respirationRecording.processSignal( ...
    processor,'respiration');

ripRecording = mpk.Recording('rip_name');
ripRecording.addSignal('rip',signal);
ripResult = ripRecording.processSignal(processor,'rip');

directResult = processor.process(signal);
localVerifyFeaturesEqual(testCase,respirationResult.Features, ...
    directResult.Features)
localVerifyFeaturesEqual(testCase,ripResult.Features,directResult.Features)
verifyEqual(testCase,respirationResult.SourceSignalLabel,'chest_belt')
verifyEqual(testCase,ripResult.SourceSignalLabel,'chest_belt')
end

function testExtractedSignalCanBeProcessedDirectly(testCase)
recording = localPhasedRecording();
processor = mpk.RespirationProcessor();
phaseSignal = recording.extractPhase('respiration','baseline');

result = processor.process(phaseSignal);

verifyClass(testCase,result,'mpk.ProcessingResult')
verifyEqual(testCase,result.SourceSignalLabel,phaseSignal.Label)
verifyEqual(testCase,result.Features(:,1), ...
    testCase.TestData.ExpectedFeatureNames)
end

function testRecordingProcessPhasesMatchesExplicitProcessing(testCase)
recording = localPhasedRecording();
processor = mpk.RespirationProcessor();
sourceSignal = recording.getSignal('respiration');
originalSamples = sourceSignal.Samples;
originalTime = sourceSignal.Time;
originalIntervals = recording.PhaseIntervalsSeconds;
originalLabels = recording.PhaseLabels;

results = recording.processPhases(processor,'respiration');

verifyClass(testCase,results,'containers.Map')
verifyEqual(testCase,sort(keys(results)),sort(originalLabels))
for index = 1:numel(originalLabels)
    phaseLabel = originalLabels{index};
    phaseSignal = recording.extractPhase('respiration',phaseLabel);
    directResult = processor.process(phaseSignal);
    mapResult = results(phaseLabel);
    verifyEqual(testCase,mapResult.Features(:,1), ...
        testCase.TestData.ExpectedFeatureNames)
    verifyFalse(testCase,any(contains(mapResult.Features(:,1),'__')))
    localVerifyFeaturesEqual(testCase,mapResult.Features, ...
        directResult.Features)
end

verifyEqual(testCase,sourceSignal.Samples,originalSamples)
verifyEqual(testCase,sourceSignal.Time,originalTime)
verifyEqual(testCase,recording.PhaseIntervalsSeconds,originalIntervals)
verifyEqual(testCase,recording.PhaseLabels,originalLabels)
end

function signal = localRespirationSignal( ...
        label,startTime,durationSeconds,fs,metadata)
if nargin < 5
    metadata = struct();
end
nSamples = round(durationSeconds * fs);
relativeTime = (0:(nSamples-1))' / fs;
time = startTime + relativeTime;
samples = sin(2*pi*0.25*relativeTime) ...
    + 0.05*sin(2*pi*0.05*relativeTime);
signal = mpk.Signal(samples,time,fs,label,'a.u.',metadata);
end

function recording = localPhasedRecording()
signal = localRespirationSignal('chest_belt',0,120,100);
recording = mpk.Recording('respiration_phases');
recording.addSignal('respiration',signal);
recording.setPhases( ...
    [0, 60; 60, 120], ...
    {'baseline','task'});
end

function localVerifyFeaturesEqual(testCase,actual,expected)
verifyEqual(testCase,actual(:,1),expected(:,1))
actualValues = cell2mat(actual(:,2));
expectedValues = cell2mat(expected(:,2));
verifyEqual(testCase,isnan(actualValues),isnan(expectedValues))
verifyEqual(testCase,isinf(actualValues),isinf(expectedValues))
finiteValues = isfinite(actualValues) & isfinite(expectedValues);
verifyEqual(testCase,actualValues(finiteValues), ...
    expectedValues(finiteValues),'AbsTol',1e-12)
verifyEqual(testCase,sign(actualValues(isinf(actualValues))), ...
    sign(expectedValues(isinf(expectedValues))))
end
