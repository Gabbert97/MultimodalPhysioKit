function tests = test_fnirs_procedural_baseline_modes
%TEST_FNIRS_PROCEDURAL_BASELINE_MODES Baseline and output extension tests.
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

function testSelfBaselineMode(testCase)
[red,infrared,time,fs,phases,ageYears] = localParentFixture();
[actualOxy,actualDeoxy,intermediates] = fNIRS_signal_processing( ...
    red,infrared,time,fs,phases,ageYears,false);
[emptyOxy,emptyDeoxy,emptyIntermediates] = fNIRS_signal_processing( ...
    red,infrared,time,fs,phases,ageYears,false, ...
    'BaselineRedValue',[],'BaselineInfraredValue',[]);

verifyEqual(testCase,emptyOxy,actualOxy)
verifyEqual(testCase,emptyDeoxy,actualDeoxy)
verifyEqual(testCase,intermediates.baseline_mode,'self')
verifyEqual(testCase,emptyIntermediates.baseline_mode,'self')
verifySize(testCase,actualOxy,[19 4])
verifySize(testCase,actualDeoxy,[19 4])
verifyEqual(testCase,actualOxy(:,1),localFeatureNames())
verifyEqual(testCase,actualDeoxy(:,1),localFeatureNames())
verifyTrue(testCase,iscell(actualOxy))
verifyTrue(testCase,iscell(actualDeoxy))
verifyTrue(testCase,all(cellfun(@isnumeric,actualOxy(:,2:end)),'all'))
verifyTrue(testCase,all(cellfun(@isnumeric,actualDeoxy(:,2:end)),'all'))

indexes = intermediates.baseline_start_index: ...
    intermediates.baseline_end_index;
verifyEqual(testCase,intermediates.baseline_red_value, ...
    mean(intermediates.filtered_red(indexes,1),1),'AbsTol',1e-14)
verifyEqual(testCase,intermediates.baseline_infrared_value, ...
    mean(intermediates.filtered_infrared(indexes,1),1),'AbsTol',1e-14)
verifyEqual(testCase,intermediates.baseline_start_index,101)
verifyEqual(testCase,intermediates.baseline_end_index,300)
end

function testPreprocessingMatchesValidatedStatements(testCase)
[red,infrared,time,fs,phases] = localParentFixture();
actual = fnirs_preprocess_currents(red,infrared,time,fs,phases);

timeFromZero = time - time(1);
[expectedInfrared,~,~] = downsampling( ...
    infrared,timeFromZero,fs,10,phases,0);
[expectedRed,expectedTime,expectedPhases] = downsampling( ...
    red,timeFromZero,fs,10,phases,0);
if expectedPhases(1,1) == 0
    expectedPhases(1,1) = 1;
end
lp = designfilt( ...
    'lowpassiir','PassbandFrequency',0.2, ...
    'StopbandFrequency',0.4,'PassbandRipple',0.5, ...
    'StopbandAttenuation',80,'SampleRate',10, ...
    'DesignMethod','butter');
expectedFilteredRed = fnirs_zero_phase_sos( ...
    lp.Coefficients,expectedRed,45);
expectedFilteredInfrared = fnirs_zero_phase_sos( ...
    lp.Coefficients,expectedInfrared,45);

verifyEqual(testCase,actual.DownsampledRed,expectedRed)
verifyEqual(testCase,actual.DownsampledInfrared,expectedInfrared)
verifyEqual(testCase,actual.Time,expectedTime)
verifyEqual(testCase,actual.PhaseIndexes,expectedPhases)
verifyEqual(testCase,actual.FilteredRed,expectedFilteredRed, ...
    'AbsTol',1e-14,'RelTol',1e-14)
verifyEqual(testCase,actual.FilteredInfrared, ...
    expectedFilteredInfrared,'AbsTol',1e-14,'RelTol',1e-14)
verifyEqual(testCase,actual.SamplingFrequency,10)
verifyEqual(testCase,actual.DownsamplingFactor,10)
end

function testExplicitZeroPhaseCurrentFilterBehavior(testCase)
lp = designfilt( ...
    'lowpassiir','PassbandFrequency',0.2, ...
    'StopbandFrequency',0.4,'PassbandRipple',0.5, ...
    'StopbandAttenuation',80,'SampleRate',10, ...
    'DesignMethod','butter');

constantSignal = 2.75*ones(600,1);
constantFiltered = fnirs_zero_phase_sos( ...
    lp.Coefficients,constantSignal,45);
verifySize(testCase,constantFiltered,size(constantSignal))
verifyTrue(testCase,all(isfinite(constantFiltered)))
verifyEqual(testCase,constantFiltered,constantSignal,'AbsTol',2e-13)

time = (0:599)'/10;
sinusoid = sin(2*pi*0.1*time);
first = fnirs_zero_phase_sos(lp.Coefficients,sinusoid,45);
second = fnirs_zero_phase_sos(lp.Coefficients,sinusoid,45);
verifyEqual(testCase,first,second)
verifyLessThan(testCase,max(abs(first)),1.1)

interior = 101:500;
[correlation,lags] = xcorr(first(interior),sinusoid(interior));
[~,maximumIndex] = max(correlation);
verifyEqual(testCase,lags(maximumIndex),0)
end

function testExplicitZeroPhaseHighpassBehavior(testCase)
hp = designfilt( ...
    'highpassiir','PassbandFrequency',0.01, ...
    'StopbandFrequency',0.006,'PassbandRipple',0.5, ...
    'StopbandAttenuation',80,'SampleRate',10, ...
    'DesignMethod','butter');
time = (0:599)'/10;
paddingLength = 598;

constant = fnirs_zero_phase_sos( ...
    hp.Coefficients,ones(size(time)),paddingLength);
trendInput = linspace(0,1,numel(time))';
trend = fnirs_zero_phase_sos( ...
    hp.Coefficients,trendInput,paddingLength);
passbandInput = sin(2*pi*0.1*time);
first = fnirs_zero_phase_sos( ...
    hp.Coefficients,passbandInput,paddingLength);
second = fnirs_zero_phase_sos( ...
    hp.Coefficients,passbandInput,paddingLength);

verifySize(testCase,constant,size(time))
verifyTrue(testCase,all(isfinite(constant)))
verifyEqual(testCase,constant,zeros(size(time)),'AbsTol',1e-13)
verifyLessThan(testCase,std(trend),0.02*std(trendInput))
verifyEqual(testCase,first,second)
verifyLessThan(testCase,max(abs(first)),1.2)

interior = 101:500;
retainedAmplitude = std(first(interior))/std(passbandInput(interior));
verifyGreaterThan(testCase,retainedAmplitude,0.98)
verifyLessThan(testCase,retainedAmplitude,1.02)
[correlation,lags] = xcorr(first(interior),passbandInput(interior));
[~,maximumIndex] = max(correlation);
verifyEqual(testCase,lags(maximumIndex),0)
end

function testDedicatedBaselineMatchesCompletePipeline(testCase)
[red,infrared,time,fs,~,ageYears] = localParentFixture();
[~,~,intermediates] = fNIRS_signal_processing( ...
    red,infrared,time,fs,[1; numel(red)],ageYears,false);
[redValue,infraredValue,diagnostics] = ...
    fnirs_compute_baseline_values(red,infrared,time,fs);

verifyEqual(testCase,redValue,intermediates.baseline_red_value, ...
    'AbsTol',1e-14)
verifyEqual(testCase,infraredValue, ...
    intermediates.baseline_infrared_value,'AbsTol',1e-14)
verifyEqual(testCase,diagnostics.Mode,'self')
verifyEqual(testCase,diagnostics.MarginSeconds,10)
verifyEqual(testCase,diagnostics.StartIndex,101)
verifyEqual(testCase,diagnostics.EndIndex,1100)
verifyEqual(testCase,diagnostics.ProcessedSamplingFrequency,10)
end

function testProvidedModeUsesValuesWithoutSelfBaseline(testCase)
[red,infrared,time,fs,~,ageYears] = localShortFixture();
phase = [1; numel(red)];
baselineRed = 10.25;
baselineInfrared = 11.5;

verifyError(testCase,@() fNIRS_signal_processing( ...
    red,infrared,time,fs,phase,ageYears,false), ...
    'MultimodalPhysioKit:ShortFnirsBaseline')
[~,~,intermediates] = fNIRS_signal_processing( ...
    red,infrared,time,fs,phase,ageYears,false, ...
    'BaselineRedValue',baselineRed, ...
    'BaselineInfraredValue',baselineInfrared);

verifyEqual(testCase,intermediates.baseline_mode,'provided')
verifyEqual(testCase,intermediates.baseline_red_value,baselineRed)
verifyEqual(testCase,intermediates.baseline_infrared_value,baselineInfrared)
verifyEmpty(testCase,intermediates.baseline_start_index)
verifyEmpty(testCase,intermediates.baseline_end_index)
end

function testBaselinePairValidation(testCase)
[red,infrared,time,fs,phases,ageYears] = localParentFixture();
baseCall = @(arguments) fNIRS_signal_processing( ...
    red,infrared,time,fs,phases,ageYears,false,arguments{:});

verifyError(testCase,@() baseCall({'BaselineRedValue',10}), ...
    'MultimodalPhysioKit:FNIRS:IncompleteBaselinePair')
verifyError(testCase,@() baseCall({'BaselineInfraredValue',11}), ...
    'MultimodalPhysioKit:FNIRS:IncompleteBaselinePair')

invalidValues = {true,[1 2],1+1i,NaN,Inf,0,-1};
for index = 1:numel(invalidValues)
    invalid = invalidValues{index};
    verifyError(testCase,@() baseCall({ ...
        'BaselineRedValue',invalid,'BaselineInfraredValue',11}), ...
        'MultimodalPhysioKit:FNIRS:InvalidBaselineRedValue')
    verifyError(testCase,@() baseCall({ ...
        'BaselineRedValue',10,'BaselineInfraredValue',invalid}), ...
        'MultimodalPhysioKit:FNIRS:InvalidBaselineInfraredValue')
end
end

function testProvidedBaselineConcentrationAndIntermediateConsistency(testCase)
[red,infrared,time,fs,phases,ageYears] = localParentFixture();
baselineRed = 10.2;
baselineInfrared = 11.3;
[oxy,deoxy,intermediates] = fNIRS_signal_processing( ...
    red,infrared,time,fs,phases,ageYears,false, ...
    'BaselineRedValue',baselineRed, ...
    'BaselineInfraredValue',baselineInfrared);
[expectedOxy,expectedDeoxy,baseline] = LambertBeer( ...
    intermediates.filtered_red,intermediates.filtered_infrared, ...
    intermediates.downsampled_time, ...
    intermediates.processed_sampling_frequency, ...
    intermediates.processed_phase_indexes,ageYears,false, ...
    'BaselineRedValue',baselineRed, ...
    'BaselineInfraredValue',baselineInfrared);
[expectedOxyFeatures,oxyDiagnostics] = features_extraction( ...
    expectedOxy,intermediates.downsampled_time, ...
    intermediates.processed_phase_indexes,10,false);
[expectedDeoxyFeatures,deoxyDiagnostics] = features_extraction( ...
    expectedDeoxy,intermediates.downsampled_time, ...
    intermediates.processed_phase_indexes,10,false);

verifyEqual(testCase,baseline.Mode,'provided')
verifyEqual(testCase,intermediates.oxy_concentration,expectedOxy, ...
    'AbsTol',1e-14,'RelTol',1e-14)
verifyEqual(testCase,intermediates.deoxy_concentration,expectedDeoxy, ...
    'AbsTol',1e-14,'RelTol',1e-14)
verifyEqual(testCase,intermediates.filtered_oxy_concentration, ...
    oxyDiagnostics.FilteredSignal,'AbsTol',1e-14,'RelTol',1e-14)
verifyEqual(testCase,intermediates.filtered_deoxy_concentration, ...
    deoxyDiagnostics.FilteredSignal,'AbsTol',1e-14,'RelTol',1e-14)
verifyEqual(testCase,oxy,expectedOxyFeatures)
verifyEqual(testCase,deoxy,expectedDeoxyFeatures)
verifyEqual(testCase,intermediates.processed_sampling_frequency,10)
verifyTrue(testCase,all(diff(intermediates.downsampled_time) > 0))
verifyEqual(testCase,numel(intermediates.filtered_red), ...
    numel(intermediates.filtered_infrared))
end

function testChangingProvidedBaselineDoesNotChangeCurrentFiltering(testCase)
[red,infrared,time,fs,phases,ageYears] = localParentFixture();
[~,~,first] = fNIRS_signal_processing( ...
    red,infrared,time,fs,phases,ageYears,false, ...
    'BaselineRedValue',10.1,'BaselineInfraredValue',11.1);
[~,~,second] = fNIRS_signal_processing( ...
    red,infrared,time,fs,phases,ageYears,false, ...
    'BaselineRedValue',10.6,'BaselineInfraredValue',11.8);

verifyEqual(testCase,first.filtered_red,second.filtered_red)
verifyEqual(testCase,first.filtered_infrared,second.filtered_infrared)
verifyNotEqual(testCase,first.oxy_concentration,second.oxy_concentration)
verifyNotEqual(testCase,first.deoxy_concentration, ...
    second.deoxy_concentration)
end

function testIndependentWindowsReuseProvidedPair(testCase)
[red,infrared,time,fs,phases,ageYears] = localParentFixture();
baselinePhase = phases(:,1);
[baselineRed,baselineInfrared] = fnirs_compute_baseline_values( ...
    red(baselinePhase(1):baselinePhase(2)), ...
    infrared(baselinePhase(1):baselinePhase(2)), ...
    time(baselinePhase(1):baselinePhase(2)),fs);

for phaseNumber = 2:3
    indexes = phases(1,phaseNumber):phases(2,phaseNumber);
    fullWindow = [1; numel(indexes)];
    [~,~,self] = fNIRS_signal_processing( ...
        red(indexes),infrared(indexes),time(indexes),fs, ...
        fullWindow,ageYears,false);
    [~,~,provided] = fNIRS_signal_processing( ...
        red(indexes),infrared(indexes),time(indexes),fs, ...
        fullWindow,ageYears,false, ...
        'BaselineRedValue',baselineRed, ...
        'BaselineInfraredValue',baselineInfrared);

    verifyEqual(testCase,provided.baseline_mode,'provided')
    verifyEqual(testCase,provided.baseline_red_value,baselineRed)
    verifyEqual(testCase,provided.baseline_infrared_value,baselineInfrared)
    verifyNotEqual(testCase,self.baseline_red_value,baselineRed)
    verifyNotEqual(testCase,self.baseline_infrared_value,baselineInfrared)
    verifyEqual(testCase,self.filtered_red,provided.filtered_red)
    verifyEqual(testCase,self.filtered_infrared,provided.filtered_infrared)
end
end

function testInputsAndNoSideEffects(testCase)
[red,infrared,time,fs,phases,ageYears] = localParentFixture();
redBefore = red;
infraredBefore = infrared;
timeBefore = time;
figuresBefore = findall(groot,'Type','figure');
directoryBefore = dir;

fNIRS_signal_processing( ...
    red,infrared,time,fs,phases,ageYears,false);

verifyEqual(testCase,red,redBefore)
verifyEqual(testCase,infrared,infraredBefore)
verifyEqual(testCase,time,timeBefore)
verifyEqual(testCase,findall(groot,'Type','figure'),figuresBefore)
directoryAfter = dir;
verifyEqual(testCase,{directoryAfter.name},{directoryBefore.name})
end

function [red,infrared,time,fs,phases,ageYears] = localParentFixture()
fs = 100;
durationSeconds = 120;
startTime = 3.7;
elapsed = (0:1/fs:durationSeconds - 1/fs)';
time = startTime + elapsed;
red = 10 + 0.08*sin(2*pi*0.025*elapsed) + ...
    0.015*sin(2*pi*0.3*elapsed) + 0.0005*elapsed;
infrared = 11 + 0.07*cos(2*pi*0.02*elapsed) + ...
    0.012*cos(2*pi*0.25*elapsed) - 0.0003*elapsed;
phaseSamples = 40*fs;
phases = [ ...
    1 phaseSamples + 1 2*phaseSamples + 1; ...
    phaseSamples 2*phaseSamples 3*phaseSamples];
ageYears = 35;
end

function [red,infrared,time,fs,phases,ageYears] = localShortFixture()
fs = 100;
durationSeconds = 20;
elapsed = (0:1/fs:durationSeconds - 1/fs)';
time = 2.4 + elapsed;
red = 10 + 0.03*sin(2*pi*0.05*elapsed);
infrared = 11 + 0.03*cos(2*pi*0.05*elapsed);
phases = [1; numel(red)];
ageYears = 35;
end

function names = localFeatureNames()
names = cellstr([ ...
    "min"; "max"; "mean"; "variance"; "std"; "skewness"; ...
    "kurtosis"; "power_band"; "max_frequency"; ...
    "median_frequency"; "spectral_entropy"; "mean_diff"; ...
    "peak_to_peak"; "slope"; "zero_crossing"; "polarity"; ...
    "entropy"; "auc"; "total_energy"]);
end
