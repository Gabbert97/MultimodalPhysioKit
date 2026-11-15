function tests = test_eda_processor
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.RepositoryRoot = repositoryRoot;
testCase.TestData.OriginalPath = path;
addpath(genpath(fullfile(repositoryRoot,'matlab')))
addpath(fullfile(repositoryRoot,'matlab','tests','helpers'))
testCase.TestData.ExpectedFeatureNames = {
    'eda_scl_mean'
    'eda_scl_std'
    'eda_scl_slope'
    'eda_scr_amplitude_mean'
    'eda_scr_amplitude_std'
    'eda_scr_rise_time_mean'
    'eda_scr_rise_time_std'
    'eda_scr_per_min'
    };
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath)
end

function testConstructorDefaultsDoNotRequireDependency(testCase)
cleanup = temporarilyHideCvxEda(); %#ok<NASGU>
processor = mpk.EDAProcessor();

verifyEqual(testCase,processor.CvxedaPath,'')
verifyEmpty(testCase,processor.DecompositionCallable)
verifyEqual(testCase,processor.CacheDirectory,'')
verifyEqual(testCase,processor.OutputLabel,'eda')
verifyEqual(testCase,processor.MinimumScrAmplitude,0.01)
verifyFalse(testCase,isprop(processor,'ShowPlots'))
verifyFalse(testCase,any(strcmp(properties(processor),'Signal')))
end

function testConstructorStoresConfiguration(testCase)
callable = @fake_cvxeda;
processor = mpk.EDAProcessor( ...
    'CvxedaPath',"external/cvxEDA", ...
    'DecompositionCallable',callable, ...
    'CacheDirectory',"cache", ...
    'OutputLabel',"skin");

verifyEqual(testCase,processor.CvxedaPath,'external/cvxEDA')
verifyEqual(testCase,processor.DecompositionCallable,callable)
verifyEqual(testCase,processor.CacheDirectory,'cache')
verifyEqual(testCase,processor.OutputLabel,'skin')
end

function testConstructorRejectsInvalidConfiguration(testCase)
verifyError(testCase,@() mpk.EDAProcessor('CvxedaPath',42), ...
    'mpk:EDAProcessor:InvalidConfiguration')
verifyError(testCase,@() mpk.EDAProcessor( ...
    'DecompositionCallable','not_callable'), ...
    'mpk:EDAProcessor:InvalidDecompositionCallable')
verifyError(testCase,@() mpk.EDAProcessor('CacheDirectory',struct()), ...
    'mpk:EDAProcessor:InvalidConfiguration')
verifyError(testCase,@() mpk.EDAProcessor('OutputLabel',{}), ...
    'mpk:EDAProcessor:InvalidConfiguration')
verifyError(testCase,@() mpk.EDAProcessor('MinimumScrAmplitude',-0.01), ...
    'mpk:EDAProcessor:InvalidMinimumScrAmplitude')
verifyError(testCase,@() mpk.EDAProcessor('MinimumScrAmplitude',NaN), ...
    'mpk:EDAProcessor:InvalidMinimumScrAmplitude')
end

function testScrTroughToPeakPairingAndThreshold(testCase)
filtered = [0 -0.005 0.005 0 -0.02 0.03 0];
[amplitudes,riseTimes,maxima,minima] = ...
    scr_pair_events(filtered,[3 6],[2 5],100);
verifyEqual(testCase,amplitudes,[0.01 0.05],'AbsTol',1e-15)
verifyEqual(testCase,riseTimes,[0.01 0.01],'AbsTol',eps)
verifyEqual(testCase,maxima,[3 6])
verifyEqual(testCase,minima,[2 5])

% A maximum without a preceding local minimum is skipped explicitly.
[amplitudes,~,maxima,minima] = ...
    scr_pair_events(filtered,[1 3],[2],100);
verifyEqual(testCase,amplitudes,0.01,'AbsTol',1e-15)
verifyEqual(testCase,maxima,3)
verifyEqual(testCase,minima,2)
end

function testScrMinimumAmplitudeAndEventPopulationPolicies(testCase)
fs = 100;
time = (0:599)' / fs;
noEvents = zeros(size(time));
[amplitudes,riseTimes,count,maxima,minima] = ...
    scr_processing(noEvents,time,fs,0.01);
verifyEmpty(testCase,amplitudes)
verifyEmpty(testCase,riseTimes)
verifyEqual(testCase,count,0)
verifyEmpty(testCase,maxima)
verifyEmpty(testCase,minima)

candidateAmplitudes = [1e-8 0.009999 0.01 0.02];
verifyEqual(testCase,scr_accept_events(candidateAmplitudes,0.01), ...
    [false false true true])
end

function testConfiguredMinimumScrAmplitudeReachesFeatureScreening(testCase)
fake_cvxeda('reset');
processor = mpk.EDAProcessor( ...
    'DecompositionCallable',@fake_cvxeda, ...
    'MinimumScrAmplitude',1.0);
features = processor.extractFeatures(localEdaSignal('eda',0,30,100));
values = cell2mat(features(:,2));
verifyEqual(testCase,values(4:8),zeros(5,1),'AbsTol',0)
end

function testDownsampleUsesExactIndexesAndTimes(testCase)
signal = localEdaSignal('eda',0,20,200);
processor = mpk.EDAProcessor();
output = processor.downsampleSignal(signal);
expectedIndexes = 1:2:signal.count();

% In-band content passes the 40 Hz anti-alias filter intact.
verifyEqual(testCase,output.Eda,signal.Samples(expectedIndexes), ...
    'AbsTol',1e-4)
verifyEqual(testCase,output.Time,signal.Time(expectedIndexes))
verifyEqual(testCase,output.SamplingFrequency,100)
end

function testDownsampleRemovesContentAboveNewNyquist(testCase)
% 100.2 Hz interference would alias to 0.2 Hz (SCR band) under plain
% decimation from 1000 Hz to 100 Hz; it must be removed.
fs = 1000;
time = (0:(60*fs-1))' / fs;
physiological = 2 + 0.3*sin(2*pi*0.05*time);
signal = mpk.Signal( ...
    physiological + 0.5*sin(2*pi*100.2*time),time,fs,'eda','uS');
output = mpk.EDAProcessor().downsampleSignal(signal);
expected = physiological(1:10:end);
edge = 5;  % exclude the zero-phase boundary transient (50 ms)
verifyEqual(testCase,output.Eda(edge+1:end-edge), ...
    expected(edge+1:end-edge),'AbsTol',1e-3)
end

function testInvalidSignalAndSamplingFrequencyAreRejected(testCase)
processor = mpk.EDAProcessor();
verifyError(testCase,@() processor.downsampleSignal([1; 2]), ...
    'mpk:EDAProcessor:InvalidSignal')
verifyError(testCase,@() processor.process(mpk.Recording('invalid')), ...
    'mpk:EDAProcessor:InvalidSignal')
verifyError(testCase,@() processor.process(struct()), ...
    'mpk:EDAProcessor:InvalidSignal')

badFrequency = localEdaSignal('eda',0,20,150);
verifyError(testCase,@() processor.downsampleSignal(badFrequency), ...
    'mpk:EDAProcessor:InvalidSamplingFrequency')
end

function testFilteringUsesDownsampledInput(testCase)
signal = localEdaSignal('eda',0,20,200);
processor = mpk.EDAProcessor();
downsampled = processor.downsampleSignal(signal);
expected = filtering(downsampled.Eda,100,0,0.5,0,5,0);
actual = processor.filterSignal(signal);

verifyEqual(testCase,actual.EdaDownsampled,downsampled.Eda)
verifyEqual(testCase,actual.Time,downsampled.Time)
verifyEqual(testCase,actual.EdaFiltered,expected,'AbsTol',1e-12)
verifyEqual(testCase,actual.SamplingFrequency,100)
end

function testTooShortFilteringHasProcessorError(testCase)
fs = 100;
time = (0:9)' / fs;
signal = mpk.Signal(1 + time,time,fs,'eda','uS');
processor = mpk.EDAProcessor();

verifyError(testCase,@() processor.filterSignal(signal), ...
    'mpk:EDAProcessor:SignalTooShortForFiltering')
end

function testCallableReceivesFilteredSignalAndDelta(testCase)
processor = localCallableProcessor();
signal = localEdaSignal('eda',0,30,200);
filtered = processor.filterSignal(signal);
decomposition = processor.decompose(signal);
fakeState = fake_cvxeda('state');

verifyEqual(testCase,fakeState.CallCount,1)
verifyEqual(testCase,fakeState.LastInput,filtered.EdaFiltered,'AbsTol',1e-12)
verifyEqual(testCase,fakeState.LastDelta,0.01,'AbsTol',eps)
verifyEqual(testCase,decomposition.SamplingFrequency,100)
verifyEqual(testCase,decomposition.EdaFiltered,filtered.EdaFiltered, ...
    'AbsTol',1e-12)
verifySize(testCase,decomposition.Scr,size(filtered.EdaFiltered))
verifySize(testCase,decomposition.Scl,size(filtered.EdaFiltered))
verifySize(testCase,decomposition.Residuals,size(filtered.EdaFiltered))
end

function testCallableHasDependencyPriority(testCase)
fake_cvxeda('reset');
processor = mpk.EDAProcessor( ...
    'DecompositionCallable',@fake_cvxeda, ...
    'CvxedaPath','path/that/does/not/exist');
processor.decompose(localEdaSignal('eda',0,20,100));
fakeState = fake_cvxeda('state');

verifyEqual(testCase,fakeState.CallCount,1)
end

function testCvxedaDirectoryAndDirectFilePathsWork(testCase)
temporaryDirectory = localTemporaryDirectory();
cleanupDirectory = onCleanup(@() localRemoveDirectory(temporaryDirectory)); %#ok<NASGU>
cvxedaFile = localWriteFakeCvxEda(temporaryDirectory);
signal = localEdaSignal('eda',0,20,100);
originalPath = path;

directoryProcessor = mpk.EDAProcessor('CvxedaPath',temporaryDirectory);
directoryResult = directoryProcessor.decompose(signal);
fileProcessor = mpk.EDAProcessor('CvxedaPath',cvxedaFile);
fileResult = fileProcessor.decompose(signal);

verifySize(testCase,directoryResult.Scr,[signal.count() 1])
verifyEqual(testCase,fileResult.Scr,directoryResult.Scr)
verifyEqual(testCase,path,originalPath)
end

function testExistingPathCvxedaWorks(testCase)
temporaryDirectory = localTemporaryDirectory();
cleanupDirectory = onCleanup(@() localRemoveDirectory(temporaryDirectory)); %#ok<NASGU>
localWriteFakeCvxEda(temporaryDirectory);
originalPath = path;
cleanupPath = onCleanup(@() path(originalPath)); %#ok<NASGU>
addpath(temporaryDirectory,'-begin')
rehash path

processor = mpk.EDAProcessor();
decomposition = processor.decompose(localEdaSignal('eda',0,20,100));

verifyNotEmpty(testCase,decomposition.Scr)
end

function testMissingDependencyOccursAfterSignalValidation(testCase)
cleanup = temporarilyHideCvxEda(); %#ok<NASGU>
processor = mpk.EDAProcessor();

verifyError(testCase,@() processor.process([1; 2; 3]), ...
    'mpk:EDAProcessor:InvalidSignal')
badFrequency = localEdaSignal('eda',0,20,150);
verifyError(testCase,@() processor.decompose(badFrequency), ...
    'mpk:EDAProcessor:InvalidSamplingFrequency')
signal = localEdaSignal('eda',0,20,100);
verifyError(testCase,@() processor.decompose(signal), ...
    'mpk:EDAProcessor:MissingCvxEDA')
end

function testProcessingOptionsAreValidated(testCase)
processor = localCallableProcessor();
signal = localEdaSignal('eda',0,20,100);

verifyError(testCase,@() processor.decompose(signal,'ShowPlots',1), ...
    'mpk:EDAProcessor:InvalidShowPlots')
verifyError(testCase,@() processor.decompose(signal,'SavePlots',1), ...
    'mpk:EDAProcessor:InvalidSavePlots')
verifyError(testCase,@() processor.decompose( ...
    signal,'OutputDirectory',42), ...
    'mpk:EDAProcessor:InvalidConfiguration')
verifyError(testCase,@() processor.decompose(signal,'OutputLabel',42), ...
    'mpk:EDAProcessor:InvalidConfiguration')
end

function testInvalidDecompositionOutputsAreRejected(testCase)
signal = localEdaSignal('eda',0,20,100);
fake_cvxeda('reset');
fake_cvxeda('short_output');
processor = mpk.EDAProcessor( ...
    'DecompositionCallable',@fake_cvxeda);
verifyError(testCase,@() processor.decompose(signal), ...
    'mpk:EDAProcessor:InvalidDecompositionOutput')

fake_cvxeda('reset');
fake_cvxeda('nonfinite_output');
processor = mpk.EDAProcessor( ...
    'DecompositionCallable',@fake_cvxeda);
verifyError(testCase,@() processor.decompose(signal), ...
    'mpk:EDAProcessor:InvalidDecompositionOutput')
end

function testEmptyCacheDirectoryProducesNoCache(testCase)
workingDirectory = localTemporaryDirectory();
cleanupDirectory = onCleanup(@() localRemoveDirectory(workingDirectory)); %#ok<NASGU>
originalDirectory = pwd;
cleanupWorking = onCleanup(@() cd(originalDirectory)); %#ok<NASGU>
cd(workingDirectory)
processor = localCallableProcessor();
processor.decompose(localEdaSignal('eda',0,20,100));

verifyEmpty(testCase,dir('*_decomposition.mat'))
end

function testCacheReuseAcrossProcessorsAndContentInvalidation(testCase)
cacheDirectory = localTemporaryDirectory();
cleanupDirectory = onCleanup(@() localRemoveDirectory(cacheDirectory)); %#ok<NASGU>
signal = localEdaSignal('eda',0,20,100);
fake_cvxeda('reset');
first = mpk.EDAProcessor( ...
    'DecompositionCallable',@fake_cvxeda, ...
    'CacheDirectory',cacheDirectory, ...
    'OutputLabel','eda');
first.decompose(signal);
first.decompose(signal);
fakeState = fake_cvxeda('state');
verifyEqual(testCase,fakeState.CallCount,1)

fake_cvxeda('reset');
second = mpk.EDAProcessor( ...
    'DecompositionCallable',@fake_cvxeda, ...
    'CacheDirectory',cacheDirectory, ...
    'OutputLabel','eda');
second.decompose(signal);
fakeState = fake_cvxeda('state');
verifyEqual(testCase,fakeState.CallCount,0)

changedSamples = mpk.Signal( ...
    signal.Samples + 0.001,signal.Time,100,'eda','uS');
second.decompose(changedSamples);
changedTime = mpk.Signal( ...
    signal.Samples,signal.Time + 0.25,100,'eda','uS');
second.decompose(changedTime);
changedFrequency = mpk.Signal( ...
    signal.Samples,signal.Time,200,'eda','uS');
second.decompose(changedFrequency);
fakeState = fake_cvxeda('state');
verifyEqual(testCase,fakeState.CallCount,3)

cacheFiles = dir(fullfile(cacheDirectory,'eda_*_decomposition.mat'));
verifyEqual(testCase,numel(cacheFiles),4)
verifyTrue(testCase,all(cellfun( ...
    @(name) ~isempty(regexp(name, ...
        '^eda_[0-9a-f]{16}_decomposition\.mat$','once')), ...
    {cacheFiles.name})))
end

function testCacheIsNotReusedByDifferentDecomposition(testCase)
cacheDirectory = localTemporaryDirectory();
cleanupDirectory = onCleanup(@() localRemoveDirectory(cacheDirectory)); %#ok<NASGU>
signal = localEdaSignal('eda',0,20,100);
fake_cvxeda('reset');
first = mpk.EDAProcessor( ...
    'DecompositionCallable',@fake_cvxeda, ...
    'CacheDirectory',cacheDirectory, ...
    'OutputLabel','eda');
first.decompose(signal);

fake_cvxeda('reset');
second = mpk.EDAProcessor( ...
    'DecompositionCallable',@(eda,delta) fake_cvxeda(eda,delta), ...
    'CacheDirectory',cacheDirectory, ...
    'OutputLabel','eda');
second.decompose(signal);
fakeState = fake_cvxeda('state');
verifyEqual(testCase,fakeState.CallCount,1)
verifyEqual(testCase, ...
    numel(dir(fullfile(cacheDirectory,'eda_*_decomposition.mat'))),2)
end

function testMalformedCacheIsRejected(testCase)
cacheDirectory = localTemporaryDirectory();
cleanupDirectory = onCleanup(@() localRemoveDirectory(cacheDirectory)); %#ok<NASGU>
processor = localCallableProcessor(cacheDirectory);
signal = localEdaSignal('eda',0,20,100);
processor.decompose(signal);
cacheFile = dir(fullfile(cacheDirectory,'*_decomposition.mat'));
bad = 1; %#ok<NASGU>
save(fullfile(cacheDirectory,cacheFile(1).name),'bad')

verifyError(testCase,@() processor.decompose(signal), ...
    'mpk:EDAProcessor:InvalidCache')
end

function testExtractFeaturesAndProcessResult(testCase)
processor = localCallableProcessor();
signal = localEdaSignal('skin_conductance',0,60,100);
features = processor.extractFeatures(signal);
result = processor.process(signal);

verifySize(testCase,features,[8 2])
verifyEqual(testCase,features(:,1),testCase.TestData.ExpectedFeatureNames)
verifyTrue(testCase,all(cellfun( ...
    @(value) isnumeric(value) && isscalar(value) && isfinite(value), ...
    features(:,2))))
verifyFalse(testCase,any(contains(features(:,1),'__')))
localVerifyFeaturesEqual(testCase,result.Features,features)
verifyClass(testCase,result,'mpk.ProcessingResult')
verifyEqual(testCase,result.Modality,'eda')
verifyEqual(testCase,result.ProcessorName,'process_eda')
verifyEqual(testCase,result.SourceSignalLabel,'skin_conductance')
verifyEqual(testCase,result.Options,struct())
verifyEqual(testCase,result.Metadata.source,'process_eda')
verifyEqual(testCase,result.Metadata.original_sampling_frequency,100)
verifyEqual(testCase,result.Metadata.processed_sampling_frequency,100)
verifyEqual(testCase,result.Metadata.original_n_samples,6000)
verifyEqual(testCase,result.Metadata.processed_n_samples,6000)
verifyEqual(testCase,result.Metadata.start_time_seconds,signal.Time(1))
verifyEqual(testCase,result.Metadata.end_time_seconds,signal.Time(end))
verifyEqual(testCase,result.Metadata.cvxeda_source,'callable')
verifyEqual(testCase,result.Metadata.minimum_scr_amplitude_us,0.01)
verifyFalse(testCase,any(contains(fieldnames(result.Metadata),'phase', ...
    'IgnoreCase',true)))
end

function testProcessorDelegatesToCanonicalProceduralPipeline(testCase)
signal = localEdaSignal('skin_conductance',4.2,60,200);
fake_cvxeda('reset');
[expectedFeatures,expectedIntermediates] = ...
    eda_signal_processing( ...
        signal.Samples, ...
        signal.Time, ...
        [1; signal.count()], ...
        signal.SamplingFrequency, ...
        'DecompositionCallable',@fake_cvxeda, ...
        'CacheDirectory','', ...
        'OutputLabel','eda', ...
        'SignalLabel',signal.Label);
proceduralState = fake_cvxeda('state');
verifyEqual(testCase,proceduralState.CallCount,1)
verifyEqual(testCase,fieldnames(expectedIntermediates),{
    'eda_downsampled'
    'eda_filtered'
    'time'
    'sampling_frequency'
    'scr'
    'scl'
    'residuals'
    'cvxeda_source'
    })

fake_cvxeda('reset');
processor = mpk.EDAProcessor( ...
    'DecompositionCallable',@fake_cvxeda, ...
    'CacheDirectory','', ...
    'OutputLabel','eda');
actual = processor.process(signal);
processorState = fake_cvxeda('state');

verifyEqual(testCase,processorState.CallCount,1)
localVerifyFeaturesEqual(testCase,actual.Features,expectedFeatures)
verifyEqual(testCase,actual.Metadata.processed_sampling_frequency, ...
    expectedIntermediates.sampling_frequency)
verifyEqual(testCase,actual.Metadata.processed_n_samples, ...
    numel(expectedIntermediates.eda_downsampled))
verifyEqual(testCase,actual.Metadata.cvxeda_source, ...
    expectedIntermediates.cvxeda_source)
end

function testEmptyLabelFallbackAndSignalImmutability(testCase)
processor = localCallableProcessor();
metadata = struct('Sensor','synthetic');
signal = localEdaSignal('',5,30,100,metadata);
originalSamples = signal.Samples;
originalTime = signal.Time;
originalMetadata = signal.Metadata;
result = processor.process(signal);

verifyEqual(testCase,result.SourceSignalLabel,'eda')
verifyEqual(testCase,signal.Samples,originalSamples)
verifyEqual(testCase,signal.Time,originalTime)
verifyEqual(testCase,signal.SamplingFrequency,100)
verifyEqual(testCase,signal.Label,'')
verifyEqual(testCase,signal.Units,'uS')
verifyEqual(testCase,signal.Metadata,originalMetadata)
end

function testProcessorReusesConfigurationAcrossSignals(testCase)
processor = localCallableProcessor();
first = processor.process(localEdaSignal('first',0,30,100));
second = processor.process(localEdaSignal('second',2,40,100));

verifyEqual(testCase,first.SourceSignalLabel,'first')
verifyEqual(testCase,second.SourceSignalLabel,'second')
verifyEqual(testCase,first.Metadata.original_n_samples,3000)
verifyEqual(testCase,second.Metadata.original_n_samples,4000)
fakeState = fake_cvxeda('state');
verifyEqual(testCase,fakeState.CallCount,2)
end

function testPlottingAndSavingDoNotChangeFeatures(testCase)

processor = localCallableProcessor();
signal = localEdaSignal('eda',0,30,100);

previousFigures = findall(groot,'Type','figure');

% Default processing must not create figures.
expected = processor.extractFeatures(signal);
currentFigures = findall(groot,'Type','figure');
localVerifySameFigures(testCase,currentFigures,previousFigures)

% Explicit plotting must create at least one figure without changing features.
visibleFeatures = processor.extractFeatures( ...
    signal, ...
    'ShowPlots',true);

currentFigures = findall(groot,'Type','figure');
newFigures = setdiff(currentFigures,previousFigures);

cleanupFigures = onCleanup(@() localCloseFigures(newFigures)); %#ok<NASGU>

verifyNotEmpty(testCase,newFigures)
localVerifyFeaturesEqual(testCase,visibleFeatures,expected)

% Saving must create the expected file without leaving additional figures.
outputDirectory = localTemporaryDirectory();
cleanupDirectory = onCleanup(@() ...
    localRemoveDirectory(outputDirectory)); %#ok<NASGU>

actual = processor.extractFeatures( ...
    signal, ...
    'ShowPlots',false, ...
    'SavePlots',true, ...
    'OutputDirectory',outputDirectory, ...
    'OutputLabel','diagnostic');

localVerifyFeaturesEqual(testCase,actual,expected)

verifyTrue(testCase,isfile(fullfile( ...
    outputDirectory, ...
    'diagnostic_eda_decomposition.png')))

figuresAfterSaving = findall(groot,'Type','figure');
expectedFigures = [previousFigures(:); newFigures(:)];

localVerifySameFigures( ...
    testCase, ...
    figuresAfterSaving, ...
    expectedFigures)

verifyError(testCase,@() ...
    processor.process(signal,'SavePlots',true), ...
    'mpk:EDAProcessor:MissingOutputDirectory')

end

function testRecordingWholeSignalAndPhaseIntegration(testCase)
processor = localCallableProcessor();
signal = localEdaSignal('skin_conductance',0,120,100);
recording = mpk.Recording('eda_recording');
recording.addSignal('eda',signal);
recording.setPhases([0, 60; 60, 120],{'segment one','segment-two'});
processArguments = {'ShowPlots',false};

wholeResult = recording.processSignal( ...
    processor,'eda','ProcessArguments',processArguments);
directWhole = processor.process(signal,'ShowPlots',false);
localVerifyFeaturesEqual(testCase,wholeResult.Features,directWhole.Features)

phaseResults = recording.processPhases( ...
    processor,'eda','ProcessArguments',processArguments);
verifyEqual(testCase,sort(keys(phaseResults)), ...
    sort(recording.PhaseLabels))
for index = 1:numel(recording.PhaseLabels)
    label = recording.PhaseLabels{index};
    extracted = recording.extractPhase('eda',label);
    direct = processor.process(extracted,'ShowPlots',false);
    localVerifyFeaturesEqual(testCase, ...
        phaseResults(label).Features,direct.Features)
end
end

function processor = localCallableProcessor(cacheDirectory)
if nargin < 1
    cacheDirectory = '';
end
fake_cvxeda('reset');
processor = mpk.EDAProcessor( ...
    'DecompositionCallable',@fake_cvxeda, ...
    'CacheDirectory',cacheDirectory);
end

function signal = localEdaSignal(label,startTime,durationSeconds,fs,metadata)
if nargin < 5
    metadata = struct();
end
nSamples = round(durationSeconds * fs);
relativeTime = (0:(nSamples-1))' / fs;
time = startTime + relativeTime;
samples = 2 + 0.02*sin(2*pi*0.05*relativeTime) ...
    + 0.002*cos(2*pi*0.2*relativeTime);
signal = mpk.Signal(samples,time,fs,label,'uS',metadata);
end

function cleanup = temporarilyHideCvxEda()
originalPath = path;
cleanup = onCleanup(@() localRestorePath(originalPath));
clear cvxEDA
cvxedaFile = which('cvxEDA');
while ~isempty(cvxedaFile)
    rmpath(fileparts(cvxedaFile))
    rehash path
    clear cvxEDA
    cvxedaFile = which('cvxEDA');
end
end

function directory = localTemporaryDirectory()
directory = tempname;
mkdir(directory)
end

function localRemoveDirectory(directory)
if isfolder(directory)
    rmdir(directory,'s')
end
end

function localRestorePath(originalPath)
path(originalPath)
rehash path
end

function filePath = localWriteFakeCvxEda(directory)
filePath = fullfile(directory,'cvxEDA.m');
fileIdentifier = fopen(filePath,'w');
if fileIdentifier < 0
    error('test:UnableToCreateFakeCvxEDA', ...
        'Unable to create fake cvxEDA.m.')
end
cleanup = onCleanup(@() fclose(fileIdentifier)); %#ok<NASGU>
fprintf(fileIdentifier,[ ...
    'function [scr,o2,scl,o4,o5,residuals,o7] = cvxEDA(eda,delta,varargin)\n' ...
    'n = numel(eda); t = (0:(n-1))'' .* delta;\n' ...
    'scr = 0.2 + 0.1*sin(2*pi*0.2*t);\n' ...
    'scl = 1.5 + 0.0002*(0:(n-1))'';\n' ...
    'residuals = zeros(n,1);\n' ...
    'o2 = zeros(n,1); o4 = o2; o5 = o2; o7 = o2;\n' ...
    'end\n']);
end

function localVerifyFeaturesEqual(testCase,actual,expected)
verifyEqual(testCase,actual(:,1),expected(:,1))
verifyEqual(testCase,cell2mat(actual(:,2)), ...
    cell2mat(expected(:,2)),'AbsTol',1e-12)
end

function localVerifySameFigures(testCase,actualFigures,expectedFigures)

verifyEmpty( ...
    testCase, ...
    setdiff(actualFigures,expectedFigures))

verifyEmpty( ...
    testCase, ...
    setdiff(expectedFigures,actualFigures))

end

function localCloseFigures(figures)

validFigures = figures(isgraphics(figures,'figure'));

if ~isempty(validFigures)
    close(validFigures)
end

end
