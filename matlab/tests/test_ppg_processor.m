function tests = test_ppg_processor
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

function testConstantRateRecoversKnownPulseShapeAndHeartRate(testCase)
signal = localPpg(repmat(0.8,1,150));          % 75 bpm, 120 s
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[]).process(signal);
verifyEqual(testCase,localFeature(result,'ppg_bpm_mean'),75,'AbsTol',1e-6)
verifyEqual(testCase,localFeature(result,'ppg_bpm_std'),0,'AbsTol',1e-6)
% The 0.5 Hz high-pass of the cleaned signal trims amplitude slightly (~0.2 %).
verifyEqual(testCase,localFeature(result,'ppg_pulse_amplitude_mean'),1,'RelTol',0.01)
verifyEqual(testCase,localFeature(result,'ppg_rise_time_mean'),0.15,'AbsTol',0.01)
verifyEqual(testCase,localFeature(result,'ppg_pulse_duration_mean'),0.8,'AbsTol',0.01)
verifyEqual(testCase,localFeature(result,'ppg_pnn50'),0)
end

function testDetectsHeartRatesAbove90Bpm(testCase)
% 150 bpm with a 1 Hz interference at 10 % of the pulse amplitude.
clean = localPpg(repmat(0.4,1,300));
signal = mpk.Signal(clean.Samples + 0.1*sin(2*pi*1.0*clean.Time), ...
    clean.Time,100,'ppg','a.u.');
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[]).process(signal);
verifyEqual(testCase,localFeature(result,'ppg_bpm_mean'),150,'AbsTol',1e-6)
end

function testPnn50OfAlternatingIntervalsIsOne(testCase)
signal = localPpg(repmat([0.8 0.9],1,70));
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[]).process(signal);
verifyEqual(testCase,localFeature(result,'ppg_pnn50'),1,'AbsTol',1e-12)
end

function testMalformedPulseIsRejectedAndItsIntervalsAreSkipped(testCase)
rises = repmat(0.15,1,100);
rises(51) = 0.6;                                % peak in the second half
signal = localPpg(repmat(0.8,1,100),rises);
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[], ...
    'ReturnIntermediates',true).process(signal);
valid = result.Intermediates.pulse_valid;
onsets = result.Intermediates.pulse_onset_indices;
rejectedOnsetTimes = signal.Time(onsets(~valid));
verifyEqual(testCase,rejectedOnsetTimes,50*0.8,'AbsTol',0.05)
verifyEqual(testCase,result.Metadata.n_accepted_ibis,numel(valid) - 1 - 2)
verifyEqual(testCase,localFeature(result,'ppg_bpm_std'),0,'AbsTol',1e-6)
end

function testPulseShapeIsMeasuredOnDriftFreeCleanedSignal(testCase)
% Slow drift 8x the pulse amplitude must not distort pulse characteristics.
clean = localPpg(repmat(0.8,1,150));
signal = mpk.Signal(clean.Samples + 8*sin(2*pi*0.05*clean.Time), ...
    clean.Time,100,'ppg','a.u.');
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[], ...
    'ReturnIntermediates',true).process(signal);
verifyGreaterThanOrEqual(testCase,result.Metadata.n_valid_pulses, ...
    result.Metadata.n_pulses - 2)
verifyEqual(testCase,localFeature(result,'ppg_rise_time_mean'),0.15,'AbsTol',0.02)
verifyEqual(testCase,localFeature(result,'ppg_pulse_amplitude_mean'),1,'RelTol',0.05)
verifyEqual(testCase,numel(result.Intermediates.cleaned_signal),signal.count())
end

function testShortRealTimeWindowWithoutFrequencyAnalysis(testCase)
signal = localPpg(repmat(0.8,1,13));            % ~10 s window
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[]).process(signal);
verifySize(testCase,result.Features,[14 2])
verifyError(testCase,@() mpk.PPGProcessor().process(signal), ...
    'MultimodalPhysioKit:PPG:TooShortForFrequencyAnalysis')
end

function testRespiratorySinusArrhythmiaIsHighFrequencyPower(testCase)
signal = localPpg(localRsaIbis(300));
for mode = [1 2]
    result = mpk.PPGProcessor('LabelFrequencyAnalysis',mode).process(signal);
    verifySize(testCase,result.Features,[20 2])
    verifyGreaterThan(testCase,localFeature(result,'ppg_phf_mean'),0.8)
    verifyLessThan(testCase,localFeature(result,'ppg_lf_hf_mean'),0.25)
end
end

function testInvalidInputsRaiseClearErrors(testCase)
verifyError(testCase,@() mpk.PPGProcessor('LabelFrequencyAnalysis',3), ...
    'mpk:PPGProcessor:InvalidFrequencyAnalysis')
verifyError(testCase,@() mpk.PPGProcessor().process(zeros(100,1)), ...
    'mpk:PPGProcessor:InvalidSignal')
slow = mpk.Signal(zeros(100,1),(0:99)'/5,5,'ppg','a.u.');
verifyError(testCase,@() mpk.PPGProcessor('LabelFrequencyAnalysis',[]).process(slow), ...
    'MultimodalPhysioKit:PPG:InvalidSamplingFrequency')
end

function value = localFeature(result,name)
value = result.Features{strcmp(result.Features(:,1),name),2};
end

function signal = localPpg(ibis,rises)
% Synthetic PPG: foot at 0, smooth rise to 1 at the rise time, decay to 0.
fs = 100;
if nargin < 2
    rises = repmat(0.15,size(ibis));
end
samples = [];
for k = 1:numel(ibis)
    n = round(ibis(k)*fs);
    r = round(rises(k)*fs);
    tau = (0:n-1)';
    pulse = (1 - (tau - r)/(n - r)).^2;
    pulse(tau < r) = 0.5 - 0.5*cos(pi*tau(tau < r)/r);
    samples = [samples; pulse]; %#ok<AGROW>
end
samples = [samples; 0];
time = (0:numel(samples)-1)'/fs;
signal = mpk.Signal(samples,time,fs,'ppg','a.u.');
end

function ibis = localRsaIbis(durationSeconds)
ibis = [];
elapsed = 0;
while elapsed < durationSeconds
    ibis(end+1) = 0.9 + 0.05*sin(2*pi*0.25*elapsed); %#ok<AGROW>
    elapsed = elapsed + ibis(end);
end
end
