function tests = test_fnirs_processor_smoke
%TEST_FNIRS_PROCESSOR_SMOKE Signal-based fNIRS processor integration tests.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.OriginalPath = path;
addpath(genpath(fullfile(repositoryRoot,'matlab')))
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath)
close all force
end

function testSelfAndProvidedProceduralParity(testCase)
[red,infrared] = localSignals(40,100,4.2);
ageYears = 35;
fullPhase = [1; red.count()];
processor = mpk.FNIRSProcessor( ...
    'AgeYears',ageYears,'ReturnIntermediates',true);

[expectedSelfOxy,expectedSelfDeoxy,expectedSelfIntermediates] = ...
    fNIRS_signal_processing( ...
        red.Samples,infrared.Samples,red.Time, ...
        red.SamplingFrequency,fullPhase,ageYears,false);
selfResult = processor.process(red,infrared);
verifyFeatureGroups(testCase,selfResult.Features, ...
    expectedSelfOxy,expectedSelfDeoxy)
verifyEqual(testCase,selfResult.Intermediates,expectedSelfIntermediates)
verifyEqual(testCase,selfResult.Metadata.baseline_mode,'self')

baselineRed = expectedSelfIntermediates.baseline_red_value;
baselineInfrared = expectedSelfIntermediates.baseline_infrared_value;
[expectedProvidedOxy,expectedProvidedDeoxy, ...
    expectedProvidedIntermediates] = fNIRS_signal_processing( ...
        red.Samples,infrared.Samples,red.Time, ...
        red.SamplingFrequency,fullPhase,ageYears,false, ...
        'BaselineRedValue',baselineRed, ...
        'BaselineInfraredValue',baselineInfrared);
providedResult = processor.process( ...
    red,infrared, ...
    'BaselineRedValue',baselineRed, ...
    'BaselineInfraredValue',baselineInfrared);
verifyFeatureGroups(testCase,providedResult.Features, ...
    expectedProvidedOxy,expectedProvidedDeoxy)
verifyEqual(testCase,providedResult.Intermediates, ...
    expectedProvidedIntermediates)
verifyEqual(testCase,providedResult.Metadata.baseline_mode,'provided')
verifyEqual(testCase,providedResult.Metadata.baseline_red_value, ...
    baselineRed)
verifyEqual(testCase,providedResult.Metadata.baseline_infrared_value, ...
    baselineInfrared)
verifyEmpty(testCase,providedResult.Intermediates.baseline_start_index)
verifyEmpty(testCase,providedResult.Intermediates.baseline_end_index)
end

function testResultContractLabelsOptionsAndMetadata(testCase)
[red,infrared] = localSignals(40,100,4.2);
processor = mpk.FNIRSProcessor( ...
    'AgeYears',31,'ReturnIntermediates',true);
result = processor.process(red,infrared);

verifyClass(testCase,result,'mpk.ProcessingResult')
verifyEqual(testCase,result.Modality,'fnirs')
verifyEqual(testCase,result.ProcessorName,'process_fnirs')
verifyEqual(testCase,cellstr(result.SourceSignalLabel), ...
    {'red_current';'infrared_current'})
verifyEqual(testCase,result.Options.AgeYears,31)
verifyTrue(testCase,result.Options.ReturnIntermediates)
verifyEmpty(testCase,result.Options.BaselineRedValue)
verifyEmpty(testCase,result.Options.BaselineInfraredValue)
verifyEqual(testCase,result.Metadata.source,'process_fnirs')
verifyEqual(testCase,result.Metadata.age_years,31)
verifyEqual(testCase,result.Metadata.original_sampling_frequency,100)
verifyEqual(testCase,result.Metadata.processed_sampling_frequency,10)
verifyEqual(testCase,result.Metadata.original_n_samples,red.count())
verifyEqual(testCase,result.Metadata.processed_n_samples, ...
    numel(result.Intermediates.downsampled_time))
verifyEqual(testCase,result.Metadata.start_time_seconds,red.Time(1))
verifyEqual(testCase,result.Metadata.end_time_seconds,red.Time(end))
verifyEqual(testCase,result.Metadata.analysis_scope,'single_signal_pair')

emptyRed = mpk.Signal( ...
    red.Samples,red.Time,red.SamplingFrequency,'','uA');
emptyInfrared = mpk.Signal( ...
    infrared.Samples,infrared.Time, ...
    infrared.SamplingFrequency,'','uA');
fallback = processor.process(emptyRed,emptyInfrared);
verifyEqual(testCase,cellstr(fallback.SourceSignalLabel), ...
    {'fnirs_red';'fnirs_infrared'})
end

function testExactFeatureRepresentation(testCase)
[red,infrared] = localSignals(40,100,0);
processor = mpk.FNIRSProcessor('AgeYears',25);
result = processor.process(red,infrared);
expectedNames = localFeatureNames();

verifyTrue(testCase,isstruct(result.Features))
verifyEqual(testCase,fieldnames(result.Features),{'Oxy';'Deoxy'})
verifySize(testCase,result.Features.Oxy,[19 2])
verifySize(testCase,result.Features.Deoxy,[19 2])
verifyEqual(testCase,result.Features.Oxy(:,1),expectedNames)
verifyEqual(testCase,result.Features.Deoxy(:,1),expectedNames)
end

function testReturnIntermediatesDoesNotChangeFeatures(testCase)
[red,infrared] = localSignals(40,100,0);
without = mpk.FNIRSProcessor( ...
    'AgeYears',25,'ReturnIntermediates',false).process( ...
    red,infrared);
with = mpk.FNIRSProcessor( ...
    'AgeYears',25,'ReturnIntermediates',true).process( ...
    red,infrared);

verifyEqual(testCase,without.Intermediates,struct())
verifyEqual(testCase,without.Features,with.Features)
verifyEqual(testCase,without.Metadata.baseline_red_value, ...
    with.Metadata.baseline_red_value)
verifyEqual(testCase,without.Metadata.baseline_infrared_value, ...
    with.Metadata.baseline_infrared_value)
verifyEqual(testCase,sort(fieldnames(with.Intermediates)), ...
    sort(localIntermediateFields()))
end

function testStageMethodsAndAliases(testCase)
[red,infrared] = localSignals(40,100,4.2);
processor = mpk.FNIRSProcessor('AgeYears',35);
fullPhase = [1; red.count()];

downsampled = processor.downsampleCurrents(red,infrared);
downsampledAlias = processor.extractDownsampledCurrents(red,infrared);
[expectedInfrared,~,~] = downsampling( ...
    infrared.Samples,red.Time-red.Time(1),100,10, ...
    fullPhase,false);
[expectedRed,expectedTime,~] = downsampling( ...
    red.Samples,red.Time-red.Time(1),100,10,fullPhase,false);
verifyEqual(testCase,downsampled,downsampledAlias)
verifyEqual(testCase,downsampled.Red,expectedRed)
verifyEqual(testCase,downsampled.Infrared,expectedInfrared)
verifyEqual(testCase,downsampled.Time,expectedTime)
verifyEqual(testCase,downsampled.SamplingFrequency,10)
verifyEqual(testCase,downsampled.DownsamplingFactor,10)

filtered = processor.filterCurrents(red,infrared);
filteredAlias = processor.extractFilteredCurrents(red,infrared);
expectedFiltered = fnirs_preprocess_currents( ...
    red.Samples,infrared.Samples,red.Time,100,fullPhase);
verifyEqual(testCase,filtered,filteredAlias)
verifyEqual(testCase,filtered.Red,expectedFiltered.FilteredRed)
verifyEqual(testCase,filtered.Infrared,expectedFiltered.FilteredInfrared)

baseline = processor.computeBaselineValues(red,infrared);
[expectedRedBaseline,expectedInfraredBaseline,diagnostics] = ...
    fnirs_compute_baseline_values( ...
        red.Samples,infrared.Samples,red.Time,100);
verifyEqual(testCase,baseline.RedValue,expectedRedBaseline)
verifyEqual(testCase,baseline.InfraredValue,expectedInfraredBaseline)
verifyEqual(testCase,baseline.Mode,'self')
verifyEqual(testCase,baseline.Diagnostics,diagnostics)
end

function testConcentrationAndFeatureStageParity(testCase)
[red,infrared] = localSignals(40,100,4.2);
processor = mpk.FNIRSProcessor('AgeYears',35);
baseline = processor.computeBaselineValues(red,infrared);
arguments = { ...
    'BaselineRedValue',baseline.RedValue, ...
    'BaselineInfraredValue',baseline.InfraredValue};

concentration = processor.computeConcentrationChanges( ...
    red,infrared,arguments{:});
concentrationAlias = processor.extractConcentrationChanges( ...
    red,infrared,arguments{:});
currentData = fnirs_preprocess_currents( ...
    red.Samples,infrared.Samples,red.Time,100,[1; red.count()]);
[expectedOxy,expectedDeoxy,expectedBaseline] = LambertBeer( ...
    currentData.FilteredRed,currentData.FilteredInfrared, ...
    currentData.Time,10,currentData.PhaseIndexes,35,false, ...
    arguments{:});
verifyEqual(testCase,concentration,concentrationAlias)
verifyEqual(testCase,concentration.Oxy,expectedOxy)
verifyEqual(testCase,concentration.Deoxy,expectedDeoxy)
verifyEqual(testCase,concentration.BaselineMode,expectedBaseline.Mode)
verifyEqual(testCase,concentration.FilteredRed,currentData.FilteredRed)
verifyEqual(testCase,concentration.FilteredInfrared, ...
    currentData.FilteredInfrared)

features = processor.extractFeatures(red,infrared,arguments{:});
[expectedOxyFeatures,expectedDeoxyFeatures] = ...
    fNIRS_signal_processing( ...
        red.Samples,infrared.Samples,red.Time,100, ...
        [1; red.count()],35,false,arguments{:});
verifyFeatureGroups(testCase,features, ...
    expectedOxyFeatures,expectedDeoxyFeatures)
end

function testProvidedModeSupportsShortWindow(testCase)
[red,infrared] = localSignals(20,100,0);
processor = mpk.FNIRSProcessor( ...
    'AgeYears',25,'ReturnIntermediates',true);
verifyError(testCase,@() processor.process(red,infrared), ...
    'MultimodalPhysioKit:ShortFnirsBaseline')

result = processor.process( ...
    red,infrared, ...
    'BaselineRedValue',10.1, ...
    'BaselineInfraredValue',11.1);
verifyEqual(testCase,result.Metadata.baseline_mode,'provided')
verifyEqual(testCase,result.Metadata.baseline_red_value,10.1)
verifyEqual(testCase,result.Metadata.baseline_infrared_value,11.1)
end

function testRecordingWholeSignalAndPhaseIntegration(testCase)
[red,infrared] = localSignals(90,100,4.2);
recording = localRecording(red,infrared);
processor = mpk.FNIRSProcessor( ...
    'AgeYears',35,'ReturnIntermediates',true);

whole = recording.processSignal( ...
    processor,{'fnirs_red','fnirs_infrared'});
directWhole = processor.process(red,infrared);
verifyEqual(testCase,whole.Features,directWhole.Features)
verifyEqual(testCase,whole.Intermediates,directWhole.Intermediates)

for resetTime = [false true]
    phaseResults = recording.processPhases( ...
        processor,{'fnirs_red','fnirs_infrared'}, ...
        'ResetTime',resetTime);
    for label = {'segment_1','segment_2','segment_3'}
        redPhase = recording.extractPhase( ...
            'fnirs_red',label{1},'ResetTime',resetTime);
        infraredPhase = recording.extractPhase( ...
            'fnirs_infrared',label{1},'ResetTime',resetTime);
        expected = processor.process(redPhase,infraredPhase);
        actual = phaseResults(label{1});
        verifyEqual(testCase,actual.Features,expected.Features)
        verifyEqual(testCase,actual.Intermediates, ...
            expected.Intermediates)
        verifyEqual(testCase,actual.Metadata.start_time_seconds, ...
            redPhase.Time(1))
        verifyEqual(testCase,actual.Metadata.end_time_seconds, ...
            redPhase.Time(end))
        verifyEqual(testCase,actual.Metadata.baseline_mode,'self')
    end
end
end

function testRecordingProvidedBaselineIntegration(testCase)
[red,infrared] = localSignals(90,100,4.2);
recording = localRecording(red,infrared);
processor = mpk.FNIRSProcessor( ...
    'AgeYears',35,'ReturnIntermediates',true);
baselineRedSignal = recording.extractPhase( ...
    'fnirs_red','segment_1','ResetTime',false);
baselineInfraredSignal = recording.extractPhase( ...
    'fnirs_infrared','segment_1','ResetTime',false);
baseline = processor.computeBaselineValues( ...
    baselineRedSignal,baselineInfraredSignal);
arguments = { ...
    'BaselineRedValue',baseline.RedValue, ...
    'BaselineInfraredValue',baseline.InfraredValue};
results = recording.processPhases( ...
    processor,{'fnirs_red','fnirs_infrared'}, ...
    'Phases',{'segment_2','segment_3'}, ...
    'ResetTime',false,'ProcessArguments',arguments);

for label = {'segment_2','segment_3'}
    redPhase = recording.extractPhase( ...
        'fnirs_red',label{1},'ResetTime',false);
    infraredPhase = recording.extractPhase( ...
        'fnirs_infrared',label{1},'ResetTime',false);
    expected = processor.process( ...
        redPhase,infraredPhase,arguments{:});
    actual = results(label{1});
    verifyEqual(testCase,actual.Features,expected.Features)
    verifyEqual(testCase,actual.Intermediates,expected.Intermediates)
    verifyEqual(testCase,actual.Metadata.baseline_mode,'provided')
    verifyEqual(testCase,actual.Metadata.baseline_red_value, ...
        baseline.RedValue)
    verifyEqual(testCase,actual.Metadata.baseline_infrared_value, ...
        baseline.InfraredValue)
    verifyEmpty(testCase,actual.Intermediates.baseline_start_index)
    verifyEmpty(testCase,actual.Intermediates.baseline_end_index)
end
end

function testSignalImmutabilityReuseAndNoHiddenBaseline(testCase)
[red,infrared] = localSignals(40,100,4.2);
processor = mpk.FNIRSProcessor( ...
    'AgeYears',35,'ReturnIntermediates',true);
redSnapshot = localSignalSnapshot(red);
infraredSnapshot = localSignalSnapshot(infrared);
provided = processor.process( ...
    red,infrared,'BaselineRedValue',10.2, ...
    'BaselineInfraredValue',11.2);
self = processor.process(red,infrared);

verifyEqual(testCase,localSignalSnapshot(red),redSnapshot)
verifyEqual(testCase,localSignalSnapshot(infrared),infraredSnapshot)
verifyEqual(testCase,provided.Metadata.baseline_mode,'provided')
verifyEqual(testCase,self.Metadata.baseline_mode,'self')
verifyNotEqual(testCase,self.Metadata.baseline_red_value,10.2)
verifyNotEqual(testCase,self.Metadata.baseline_infrared_value,11.2)
verifyEqual(testCase,properties(processor), ...
    {'AgeYears';'ReturnIntermediates'})
end

function testNoInternalSideEffectsAndExternalOutputs(testCase)
[red,infrared] = localSignals(40,100,0);
processor = mpk.FNIRSProcessor( ...
    'AgeYears',25,'ReturnIntermediates',true);
figuresBefore = findall(groot,'Type','figure');
directoryBefore = dir;
result = processor.process(red,infrared);
verifyEqual(testCase,findall(groot,'Type','figure'),figuresBefore)
directoryAfter = dir;
verifyEqual(testCase,{directoryAfter.name},{directoryBefore.name})

figureHandle = figure('Visible','off');
cleanupFigure = onCleanup(@() closeIfValid(figureHandle)); %#ok<NASGU>
axesHandle = axes(figureHandle);
plot(axesHandle,result.Intermediates.downsampled_time, ...
    result.Intermediates.oxy_concentration)
verifyTrue(testCase,isgraphics(axesHandle,'axes'))

temporaryDirectory = tempname;
mkdir(temporaryDirectory)
cleanupDirectory = onCleanup( ...
    @() removeDirectoryIfPresent(temporaryDirectory)); %#ok<NASGU>
intermediates = result.Intermediates; %#ok<NASGU>
outputFile = fullfile(temporaryDirectory,'fnirs_window.mat');
save(outputFile,'intermediates')
verifyTrue(testCase,isfile(outputFile))
end

function recording = localRecording(red,infrared)
recording = mpk.Recording('fnirs_processor_test');
recording.addSignal('fnirs_red',red);
recording.addSignal('fnirs_infrared',infrared);
startTime = red.Time(1);
recording.setPhases( ...
    [startTime startTime + 30 startTime + 60; ...
     startTime + 30 startTime + 60 startTime + 90], ...
    {'segment_1','segment_2','segment_3'});
end

function [redSignal,infraredSignal] = ...
        localSignals(durationSeconds,fs,startTime)
elapsed = (0:1/fs:durationSeconds - 1/fs)';
time = startTime + elapsed;
red = 10 + 0.08*sin(2*pi*0.025*elapsed) + ...
    0.015*sin(2*pi*0.3*elapsed) + 0.0005*elapsed;
infrared = 11 + 0.07*cos(2*pi*0.02*elapsed) + ...
    0.012*cos(2*pi*0.25*elapsed) - 0.0003*elapsed;
redSignal = mpk.Signal(red,time,fs,'red_current','uA', ...
    struct('channel',1));
infraredSignal = mpk.Signal( ...
    infrared,time,fs,'infrared_current','uA', ...
    struct('channel',2));
end

function verifyFeatureGroups(testCase,actual,expectedOxy,expectedDeoxy)
verifyEqual(testCase,actual.Oxy,expectedOxy)
verifyEqual(testCase,actual.Deoxy,expectedDeoxy)
end

function names = localFeatureNames()
names = cellstr([ ...
    "min"; "max"; "mean"; "variance"; "std"; "skewness"; ...
    "kurtosis"; "power_band"; "max_frequency"; ...
    "median_frequency"; "spectral_entropy"; "mean_diff"; ...
    "peak_to_peak"; "slope"; "zero_crossing"; "polarity"; ...
    "entropy"; "auc"; "total_energy"]);
end

function fields = localIntermediateFields()
fields = { ...
    'downsampled_red'; ...
    'downsampled_infrared'; ...
    'downsampled_time'; ...
    'processed_phase_indexes'; ...
    'processed_sampling_frequency'; ...
    'filtered_red'; ...
    'filtered_infrared'; ...
    'baseline_mode'; ...
    'baseline_red_value'; ...
    'baseline_infrared_value'; ...
    'baseline_start_index'; ...
    'baseline_end_index'; ...
    'oxy_concentration'; ...
    'deoxy_concentration'; ...
    'filtered_oxy_concentration'; ...
    'filtered_deoxy_concentration'};
end

function snapshot = localSignalSnapshot(signal)
snapshot = struct( ...
    'Samples',signal.Samples, ...
    'Time',signal.Time, ...
    'SamplingFrequency',signal.SamplingFrequency, ...
    'Label',signal.Label, ...
    'Units',signal.Units, ...
    'Metadata',signal.Metadata);
end

function closeIfValid(handle)
if isgraphics(handle,'figure')
    close(handle)
end
end

function removeDirectoryIfPresent(directory)
if isfolder(directory)
    rmdir(directory,'s')
end
end
