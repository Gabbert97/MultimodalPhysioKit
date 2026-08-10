function tests = test_ecg_processor_smoke
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

function testConstructorDefaultsAndValidation(testCase)
processor = mpk.ECGProcessor();
verifyEqual(testCase,processor.LabelFrequencyAnalysis,1)
verifyEqual(testCase,processor.ReturnIntermediates,false)

methodOne = mpk.ECGProcessor('LabelFrequencyAnalysis',1);
methodTwo = mpk.ECGProcessor('LabelFrequencyAnalysis',2, ...
    'ReturnIntermediates',true);
verifyEqual(testCase,methodOne.LabelFrequencyAnalysis,1)
verifyEqual(testCase,methodTwo.LabelFrequencyAnalysis,2)
verifyEqual(testCase,methodTwo.ReturnIntermediates,true)

invalidMethods = {0,3,NaN,Inf,1 + 1i,[1 2],'1'};
for index = 1:numel(invalidMethods)
    verifyError(testCase,@() mpk.ECGProcessor( ...
        'LabelFrequencyAnalysis',invalidMethods{index}), ...
        'mpk:ECGProcessor:InvalidFrequencyAnalysis')
end
invalidFlags = {0,1,[],[true false],'true'};
for index = 1:numel(invalidFlags)
    verifyError(testCase,@() mpk.ECGProcessor( ...
        'ReturnIntermediates',invalidFlags{index}), ...
        'mpk:ECGProcessor:InvalidReturnIntermediates')
end
end

function testRejectsAnythingExceptExactlyOneSignal(testCase)
processor = mpk.ECGProcessor();
signal = localEcgSignal(65,0,'ecg');
recording = mpk.Recording('invalid-input');
recording.addSignal('ecg',signal);

verifyError(testCase,@() processor.process([1; 2]), ...
    'mpk:ECGProcessor:InvalidSignal')
verifyError(testCase,@() processor.process({signal}), ...
    'mpk:ECGProcessor:InvalidSignal')
verifyError(testCase,@() processor.process(recording), ...
    'mpk:ECGProcessor:InvalidSignal')
verifyError(testCase,@() processor.process(struct()), ...
    'mpk:ECGProcessor:InvalidSignal')
localVerifyCallFails(testCase,@() processor.process(signal,signal))
localVerifyCallFails(testCase,@() processor.process( ...
    signal,'ShowPlots',false))
end

function testProceduralParityForBothMethods(testCase)
signal = localEcgSignal(130,4.2,'ecg');
expectedNames = localExpectedFeatureNames();

for method = [1 2]
    [expectedFeatures,expectedIntermediates] = localProcedural( ...
        signal,method);
    processor = mpk.ECGProcessor( ...
        'LabelFrequencyAnalysis',method, ...
        'ReturnIntermediates',true);
    actual = processor.process(signal);

    verifyClass(testCase,actual,'mpk.ProcessingResult')
    verifySize(testCase,actual.Features,[11 2])
    verifyEqual(testCase,actual.Features(:,1),expectedNames)
    localVerifyFeatureParity(testCase,actual.Features,expectedFeatures)
    localVerifyIntermediateParity( ...
        testCase,actual.Intermediates,expectedIntermediates,signal,method)
end
end

function testResultIdentityOptionsMetadataAndLabelFallback(testCase)
signal = localEcgSignal(65,7.5,'');
processor = mpk.ECGProcessor( ...
    'LabelFrequencyAnalysis',2, ...
    'ReturnIntermediates',true);
result = processor.process(signal);

verifyEqual(testCase,result.Modality,'ecg')
verifyEqual(testCase,result.ProcessorName,'process_ecg')
verifyEqual(testCase,result.SourceSignalLabel,'ecg')
verifyEqual(testCase,result.Options,struct( ...
    'LabelFrequencyAnalysis',2, ...
    'ReturnIntermediates',true))
verifyEqual(testCase,result.Metadata.source,'process_ecg')
verifyEqual(testCase,result.Metadata.sampling_frequency, ...
    signal.SamplingFrequency)
verifyEqual(testCase,result.Metadata.n_samples,signal.count())
verifyEqual(testCase,result.Metadata.start_time_seconds,signal.Time(1))
verifyEqual(testCase,result.Metadata.end_time_seconds,signal.Time(end))
verifyEqual(testCase,result.Metadata.frequency_analysis_method,2)
verifyEqual(testCase,result.Metadata.return_intermediates,true)
verifyEqual(testCase,result.Metadata.analysis_scope,'single_signal')
verifyFalse(testCase,any(contains(fieldnames(result.Metadata), ...
    {'phase','recording'}),'all'))
end

function testReturnIntermediatesChangesOnlyDiagnostics(testCase)
signal = localEcgSignal(130,4.2,'ecg');
compactProcessor = mpk.ECGProcessor('ReturnIntermediates',false);
diagnosticProcessor = mpk.ECGProcessor('ReturnIntermediates',true);

compact = compactProcessor.process(signal);
diagnostic = diagnosticProcessor.process(signal);

verifyEqual(testCase,compact.Intermediates,struct())
verifyTrue(testCase,all(ismember( ...
    localExpectedCommonIntermediateFields(), ...
    fieldnames(diagnostic.Intermediates))))
verifyTrue(testCase,all(ismember( ...
    localExpectedStftAliasFields(), ...
    fieldnames(diagnostic.Intermediates))))
localVerifyFeatureParity(testCase,compact.Features,diagnostic.Features)
end

function testStftAliasesAreConsistentAndExternallyUsable(testCase)
signal = localEcgSignal(130,4.2,'ecg');
processor = mpk.ECGProcessor('ReturnIntermediates',true);
result = processor.process(signal);
values = result.Intermediates;
frequency = values.frequency;

verifyEqual(testCase,values.stft_matrix,frequency.Matrix)
verifyEqual(testCase,values.stft_power,frequency.Power)
verifyEqual(testCase,values.stft_frequencies,frequency.Frequencies)
verifyEqual(testCase,values.stft_times,frequency.Time)
verifyEqual(testCase,values.stft_plf,frequency.PLF)
verifyEqual(testCase,values.stft_phf,frequency.PHF)
verifyEqual(testCase,values.stft_lf_hf,frequency.LFHF)
verifyEqual(testCase,values.stft_power,abs(values.stft_matrix).^2, ...
    'AbsTol',1e-14,'RelTol',1e-14)
verifyEqual(testCase,size(values.stft_matrix,1), ...
    numel(values.stft_frequencies))
verifyEqual(testCase,size(values.stft_matrix,2), ...
    numel(values.stft_times))
verifyEqual(testCase,size(values.stft_power),size(values.stft_matrix))
verifyEqual(testCase,numel(values.stft_plf),numel(values.stft_times))
verifyEqual(testCase,numel(values.stft_phf),numel(values.stft_times))
verifyEqual(testCase,numel(values.stft_lf_hf), ...
    numel(values.stft_times))
verifyTrue(testCase,all(diff(values.stft_frequencies) > 0))
verifyTrue(testCase,all(diff(values.stft_times) > 0))

figureHandle = figure('Visible','off');
figureCleanup = onCleanup(@() localCloseFigure(figureHandle));
axesHandle = axes(figureHandle);
imagesc(axesHandle,values.stft_times,values.stft_frequencies, ...
    values.stft_power);
axis(axesHandle,'xy')
verifyTrue(testCase,isgraphics(axesHandle,'axes'))
clear figureCleanup

temporaryDirectory = tempname;
mkdir(temporaryDirectory)
directoryCleanup = onCleanup(@() rmdir(temporaryDirectory,'s'));
frequencyFile = fullfile(temporaryDirectory,'ecg_stft_window.mat');
save(frequencyFile,'frequency')
loaded = load(frequencyFile,'frequency');
verifyEqual(testCase,loaded.frequency,frequency)
clear directoryCleanup
end

function testPsdAliasesReferenceProceduralData(testCase)
signal = localEcgSignal(65,2.5,'ecg');
[~,expected] = localProcedural(signal,2);
processor = mpk.ECGProcessor( ...
    'LabelFrequencyAnalysis',2, ...
    'ReturnIntermediates',true);
result = processor.process(signal);

verifyEqual(testCase,result.Intermediates.psd, ...
    expected.frequency.PSD{1})
verifyEqual(testCase,result.Intermediates.psd_frequencies, ...
    expected.frequency.Frequencies{1})
verifyFalse(testCase,isfield(result.Intermediates,'stft_matrix'))
frequency = result.Intermediates.frequency;
verifyEqual(testCase,frequency.PLF(1),result.Features{6,2}, ...
    'AbsTol',1e-14)
verifyEqual(testCase,frequency.PHF(1),result.Features{8,2}, ...
    'AbsTol',1e-14)
verifyEqual(testCase,frequency.LFHF(1),result.Features{10,2}, ...
    'AbsTol',1e-14)
verifyEqual(testCase,result.Features{7,2},0)
verifyEqual(testCase,result.Features{9,2},0)
verifyEqual(testCase,result.Features{11,2},0)
end

function testPeakAndIntervalIntermediatesAreConsistentAndRepeatable(testCase)
signal = localEcgSignal(65,2.5,'ecg');
processor = mpk.ECGProcessor( ...
    'LabelFrequencyAnalysis',2, ...
    'ReturnIntermediates',true);
firstResult = processor.process(signal);
secondResult = processor.process(signal);
first = firstResult.Intermediates;
second = secondResult.Intermediates;

verifyEqual(testCase,first.filtered_signal,second.filtered_signal, ...
    'AbsTol',1e-14)
verifyEqual(testCase,first.candidate_r_peak_indices, ...
    second.candidate_r_peak_indices)
verifyEqual(testCase,first.candidate_r_peak_amplitudes, ...
    second.candidate_r_peak_amplitudes,'AbsTol',1e-14)
verifyTrue(testCase,all(diff(first.candidate_r_peak_indices) > 0))
verifyTrue(testCase,all(diff(first.accepted_r_peak_indices) > 0))
verifyTrue(testCase,all(ismember( ...
    first.accepted_r_peak_indices,first.candidate_r_peak_indices)))
verifyEqual(testCase,first.candidate_r_peak_amplitudes, ...
    first.filtered_signal(first.candidate_r_peak_indices), ...
    'AbsTol',1e-14)
verifyEqual(testCase,first.accepted_r_peak_amplitudes, ...
    first.filtered_signal(first.accepted_r_peak_indices), ...
    'AbsTol',1e-14)
verifyTrue(testCase,all(isfinite(first.ibi)))
verifyTrue(testCase,all(first.ibi > 0))
verifyTrue(testCase,all(isfinite(first.bpm)))
verifyEqual(testCase,numel(first.ibi),numel(first.ibi_times))
verifyEqual(testCase,numel(first.bpm),numel(first.bpm_times))
verifyEqual(testCase,first.ibi_times,first.bpm_times)
end

function testSignalImmutabilityAndProcessorReuse(testCase)
firstSignal = localEcgSignal(65,1.5,'first_ecg');
secondSignal = localEcgSignal(70,100,'second_ecg');
firstSnapshot = localSnapshot(firstSignal);
secondSnapshot = localSnapshot(secondSignal);
processor = mpk.ECGProcessor('ReturnIntermediates',true);

firstResult = processor.process(firstSignal);
secondResult = processor.process(secondSignal);

localVerifySignalUnchanged(testCase,firstSignal,firstSnapshot)
localVerifySignalUnchanged(testCase,secondSignal,secondSnapshot)
verifyEqual(testCase,firstResult.SourceSignalLabel,'first_ecg')
verifyEqual(testCase,secondResult.SourceSignalLabel,'second_ecg')
verifyEqual(testCase,firstResult.Metadata.start_time_seconds,1.5)
verifyEqual(testCase,secondResult.Metadata.start_time_seconds,100)
verifyNotEqual(testCase,firstResult.Intermediates.stft_times, ...
    secondResult.Intermediates.stft_times)
processorProperties = properties(processor);

forbiddenProperties = {
    'Signal'
    'Recording'
    'PhaseIntervals'
    };

verifyFalse( ...
    testCase, ...
    any(ismember(processorProperties,forbiddenProperties)));
end

function testRecordingProcessSignalUsesCompleteStoredSignal(testCase)
signal = localEcgSignal(65,3.0,'chest_ecg');
recording = mpk.Recording('process-signal');
recording.addSignal('channel_1',signal);
processor = mpk.ECGProcessor( ...
    'LabelFrequencyAnalysis',1, ...
    'ReturnIntermediates',true);

throughRecording = recording.processSignal(processor,'channel_1');
direct = processor.process(signal);

localVerifyFeatureParity(testCase,throughRecording.Features,direct.Features)
verifyEqual(testCase,throughRecording.Intermediates.stft_matrix, ...
    direct.Intermediates.stft_matrix)
verifyEqual(testCase,throughRecording.SourceSignalLabel,'chest_ecg')
end

function testRecordingPhasesAreIndependentForBothTimeModes(testCase)
signal = localEcgSignal(130,4.2,'ecg');
recording = mpk.Recording('phase-processing');
recording.addSignal('ecg',signal);
phaseBoundary = signal.Time(1) + 65;
signalStop = signal.Time(end) + 1/signal.SamplingFrequency;
recording.setPhases( ...
    [signal.Time(1) phaseBoundary; phaseBoundary signalStop], ...
    {'segment one','segment-two'});
processor = mpk.ECGProcessor('ReturnIntermediates',true);

absoluteResults = recording.processPhases( ...
    processor,'ecg','ResetTime',false);
resetResults = recording.processPhases( ...
    processor,'ecg','ResetTime',true);

labels = recording.PhaseLabels;
for index = 1:numel(labels)
    label = labels{index};
    absoluteSignal = recording.extractPhase( ...
        'ecg',label,'ResetTime',false);
    resetSignal = recording.extractPhase( ...
        'ecg',label,'ResetTime',true);
    absoluteDirect = processor.process(absoluteSignal);
    resetDirect = processor.process(resetSignal);
    absoluteResult = absoluteResults(label);
    resetResult = resetResults(label);

    verifyClass(testCase,absoluteResult,'mpk.ProcessingResult')
    localVerifyFeatureParity(testCase,absoluteResult.Features, ...
        absoluteDirect.Features)
    localVerifyFeatureParity(testCase,resetResult.Features, ...
        resetDirect.Features)
    verifyEqual(testCase,absoluteResult.Intermediates.stft_matrix, ...
        absoluteDirect.Intermediates.stft_matrix)
    verifyEqual(testCase,resetResult.Intermediates.stft_matrix, ...
        resetDirect.Intermediates.stft_matrix)
    verifyGreaterThanOrEqual(testCase, ...
        min(absoluteResult.Intermediates.stft_times), ...
        absoluteSignal.Time(1))
    verifyLessThanOrEqual(testCase, ...
        max(absoluteResult.Intermediates.stft_times), ...
        absoluteSignal.Time(end))
    verifyGreaterThanOrEqual(testCase, ...
        min(absoluteResult.Intermediates.frequency.UniformTime), ...
        absoluteSignal.Time(1))
    verifyLessThanOrEqual(testCase, ...
        max(absoluteResult.Intermediates.frequency.UniformTime), ...
        absoluteSignal.Time(end))
    verifyEqual(testCase,absoluteResult.Metadata.start_time_seconds, ...
        absoluteSignal.Time(1))
    verifyEqual(testCase,absoluteResult.Metadata.end_time_seconds, ...
        absoluteSignal.Time(end))
    verifyEqual(testCase,resetResult.Metadata.start_time_seconds,0)
    verifyGreaterThanOrEqual(testCase, ...
        min(resetResult.Intermediates.stft_times),0)
    verifyLessThanOrEqual(testCase, ...
        max(resetResult.Intermediates.stft_times),resetSignal.Time(end))
end

firstAbsolute = absoluteResults(labels{1});
secondAbsolute = absoluteResults(labels{2});
verifyLessThan(testCase,max(firstAbsolute.Intermediates.stft_times), ...
    min(secondAbsolute.Intermediates.stft_times))
verifyNotEqual(testCase,firstAbsolute.Intermediates.stft_matrix, ...
    secondAbsolute.Intermediates.stft_matrix)
end

function testProcessingCreatesNoFiguresOrFiles(testCase)
signal = localEcgSignal(65,0,'ecg');
processor = mpk.ECGProcessor('ReturnIntermediates',true);
temporaryDirectory = tempname;
mkdir(temporaryDirectory)
originalDirectory = pwd;
cleanup = onCleanup(@() localRestoreDirectory( ...
    originalDirectory,temporaryDirectory));
cd(temporaryDirectory)
filesBefore = dir;
figuresBefore = findall(0,'Type','figure');

processor.process(signal);

verifyEqual(testCase,findall(0,'Type','figure'),figuresBefore)
filesAfter = dir;
verifyEqual(testCase,{filesAfter.name},{filesBefore.name})
clear cleanup
end

function [features,intermediates] = localProcedural(signal,method)
[features,intermediates] = ecg_signal_processing( ...
    signal.Samples, ...
    signal.Time, ...
    [1; signal.count()], ...
    1, ...
    signal.SamplingFrequency, ...
    method, ...
    false, ...
    0, ...
    'ecg', ...
    false);
end

function signal = localEcgSignal(durationSeconds,startTime,label)
fs = 100;
elapsed = (0:1/fs:durationSeconds - 1/fs)';
time = startTime + elapsed;
phase = mod(1.2*elapsed + 0.08*sin(2*pi*0.03*elapsed),1);
samples = 0.02*sin(2*pi*1.2*elapsed) + ...
    exp(-((phase - 0.08)/0.015).^2);
signal = mpk.Signal(samples,time,fs,label,'mV', ...
    struct('Fixture','deterministic_ecg'));
end

function names = localExpectedFeatureNames()
names = { ...
    'ecg_bpm_mean'; ...
    'ecg_bpm_std'; ...
    'ecg_bpm_median'; ...
    'ecg_pnn50'; ...
    'ecg_pnn50_max'; ...
    'ecg_plf_mean'; ...
    'ecg_plf_std'; ...
    'ecg_phf_mean'; ...
    'ecg_phf_std'; ...
    'ecg_lf_hf_mean'; ...
    'ecg_lf_hf_std'};
end

function fields = localExpectedCommonIntermediateFields()
fields = { ...
    'filtered_signal'; ...
    'candidate_r_peak_indices'; ...
    'candidate_r_peak_amplitudes'; ...
    'candidate_r_peak_times'; ...
    'accepted_r_peak_indices'; ...
    'accepted_r_peak_amplitudes'; ...
    'accepted_r_peak_times'; ...
    'ibi'; ...
    'ibi_times'; ...
    'bpm'; ...
    'bpm_times'; ...
    'frequency'};
end

function fields = localExpectedStftAliasFields()
fields = { ...
    'stft_matrix'; ...
    'stft_power'; ...
    'stft_frequencies'; ...
    'stft_times'; ...
    'stft_plf'; ...
    'stft_phf'; ...
    'stft_lf_hf'};
end

function localVerifyIntermediateParity(testCase,actual,expected,signal,method)
verifyEqual(testCase,actual.filtered_signal, ...
    expected.filtered_signal{1}(:),'AbsTol',1e-14)
verifyEqual(testCase,actual.candidate_r_peak_indices, ...
    expected.candidate_r_peak_indices{1}(:))
verifyEqual(testCase,actual.candidate_r_peak_amplitudes, ...
    expected.candidate_r_peak_amplitudes{1}(:),'AbsTol',1e-14)
verifyEqual(testCase,actual.candidate_r_peak_times, ...
    signal.Time(expected.candidate_r_peak_indices{1}(:)))
verifyEqual(testCase,actual.accepted_r_peak_indices, ...
    expected.accepted_r_peak_indices{1}(:))
verifyEqual(testCase,actual.accepted_r_peak_amplitudes, ...
    expected.accepted_r_peak_amplitudes{1}(:),'AbsTol',1e-14)
verifyEqual(testCase,actual.accepted_r_peak_times, ...
    signal.Time(expected.accepted_r_peak_indices{1}(:)))
verifyEqual(testCase,actual.ibi,expected.ibi{1}(:),'AbsTol',1e-14)
verifyEqual(testCase,actual.ibi_times,expected.ibi_times{1}(:))
verifyEqual(testCase,actual.bpm,expected.bpm{1}(:),'AbsTol',1e-14)
verifyEqual(testCase,actual.bpm_times,expected.bpm_times{1}(:))
verifyEqual(testCase,actual.frequency,expected.frequency)
if method == 1
    verifyEqual(testCase,actual.stft_matrix,expected.frequency.Matrix)
    verifyEqual(testCase,actual.stft_power,expected.frequency.Power)
    verifyEqual(testCase,actual.stft_frequencies, ...
        expected.frequency.Frequencies)
    verifyEqual(testCase,actual.stft_times,expected.frequency.Time)
else
    verifyEqual(testCase,actual.psd,expected.frequency.PSD{1})
    verifyEqual(testCase,actual.psd_frequencies, ...
        expected.frequency.Frequencies{1})
end
end

function localVerifyFeatureParity(testCase,actual,expected)
verifyEqual(testCase,size(actual),size(expected))
verifyEqual(testCase,actual(:,1),expected(:,1))
actualValues = cell2mat(actual(:,2));
expectedValues = cell2mat(expected(:,2));
verifyEqual(testCase,isnan(actualValues),isnan(expectedValues))
finite = isfinite(actualValues) & isfinite(expectedValues);
verifyEqual(testCase,actualValues(finite),expectedValues(finite), ...
    'AbsTol',1e-12,'RelTol',1e-12)
end

function snapshot = localSnapshot(signal)
snapshot = struct( ...
    'Samples',signal.Samples, ...
    'Time',signal.Time, ...
    'SamplingFrequency',signal.SamplingFrequency, ...
    'Label',signal.Label, ...
    'Units',signal.Units, ...
    'Metadata',signal.Metadata);
end

function localVerifySignalUnchanged(testCase,signal,snapshot)
verifyEqual(testCase,signal.Samples,snapshot.Samples)
verifyEqual(testCase,signal.Time,snapshot.Time)
verifyEqual(testCase,signal.SamplingFrequency,snapshot.SamplingFrequency)
verifyEqual(testCase,signal.Label,snapshot.Label)
verifyEqual(testCase,signal.Units,snapshot.Units)
verifyEqual(testCase,signal.Metadata,snapshot.Metadata)
end

function localCloseFigure(figureHandle)
if isgraphics(figureHandle)
    close(figureHandle)
end
end

function localRestoreDirectory(originalDirectory,temporaryDirectory)
cd(originalDirectory)
if isfolder(temporaryDirectory)
    rmdir(temporaryDirectory,'s')
end
end

function localVerifyCallFails(testCase,callable)
didFail = false;
try
    callable();
catch errorCaught
    didFail = true;
    verifyNotEmpty(testCase,errorCaught.identifier)
end
verifyTrue(testCase,didFail)
end
