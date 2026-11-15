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
% The 0.5-8 Hz cleaning trims amplitude slightly (~1 %) and rounds the sharp
% synthetic foot, which moves it up to 2 samples earlier.
verifyEqual(testCase,localFeature(result,'ppg_pulse_amplitude_mean'),1,'RelTol',0.02)
verifyEqual(testCase,localFeature(result,'ppg_rise_time_mean'),0.15,'AbsTol',0.02)
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
% Neighbouring beats may shift by one sample (band-pass spreads the bad pulse).
verifyLessThan(testCase,localFeature(result,'ppg_bpm_std'),0.5)
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

function testIbisAreMeasuredBetweenSystolicPeaks(testCase)
% Feet every 0.8 s, rise time alternating 0.15 / 0.25 s: systolic peaks are
% 0.9 / 0.7 s apart, and IBIs must be these peak-to-peak intervals.
signal = localPpg(repmat(0.8,1,150),repmat([0.15 0.25],1,75));
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[], ...
    'ReturnIntermediates',true).process(signal);
peakTimes = signal.Time(result.Intermediates.pulse_peak_indices);
expected = diff(peakTimes);
ibi = result.Intermediates.ibi;
verifyEqual(testCase,ibi,expected(1:numel(ibi)),'AbsTol',1e-9)
% The 8 Hz low-pass can move these asymmetric peaks by one sample (10 ms).
verifyEqual(testCase,ibi(1:2:end),repmat(ibi(1),ceil(numel(ibi)/2),1),'AbsTol',1e-9)
verifyEqual(testCase,sort(ibi(1:2)),[0.7; 0.9],'AbsTol',0.01 + 1e-9)
verifyEqual(testCase,localFeature(result,'ppg_pnn50'),1,'AbsTol',1e-12)
end

function testMaxPeakPositionControlsLatePeakRule(testCase)
% Rounded, well-formed pulses with the peak at 55 % of the pulse.
fs = 100; n = 80; r = 44;
tau = (0:n-1)';
pulse = 0.5 + 0.5*cos(pi*(tau - r)/(n - r));
pulse(tau < r) = 0.5 - 0.5*cos(pi*tau(tau < r)/r);
samples = [repmat(pulse,150,1); 0];
signal = mpk.Signal(samples,(0:numel(samples)-1)'/fs,fs,'ppg','a.u.');
verifyError(testCase,@() mpk.PPGProcessor('LabelFrequencyAnalysis',[]).process(signal), ...
    'MultimodalPhysioKit:PPG:TooFewPulses')
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[], ...
    'MaxPeakPosition',0.65).process(signal);
verifyEqual(testCase,result.Metadata.n_valid_pulses,result.Metadata.n_pulses)
% Rounded peaks can land one sample apart, so IBIs are exact to one sample.
verifyEqual(testCase,localFeature(result,'ppg_bpm_mean'),75,'AbsTol',0.1)
verifyError(testCase,@() mpk.PPGProcessor('MaxPeakPosition',1.5), ...
    'mpk:PPGProcessor:InvalidMaxPeakPosition')
end

function testDicroticWaveDoesNotSplitThePulse(testCase)
% Systolic wave at 0.15 s and a dicrotic wave half as tall at 0.45 s, one
% beat every 0.8 s: the notch between them is not a pulse foot.
signal = localTwoHumpPpg(1.0,0.5,0.15,0.45);
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[], ...
    'ReturnIntermediates',true).process(signal);
verifyGreaterThanOrEqual(testCase,result.Metadata.n_pulses,145)
verifyGreaterThanOrEqual(testCase,result.Metadata.n_valid_pulses, ...
    result.Metadata.n_pulses - 2)
verifyEqual(testCase,result.Intermediates.ibi, ...
    repmat(0.8,size(result.Intermediates.ibi)),'AbsTol',0.011)
end

function testSystolicPeakIsTheFirstHumpEvenWhenALaterOneIsTaller(testCase)
% Wrist-like pulse: systolic hump (0.85) at 0.15 s, taller late hump (1.0)
% at 0.32 s. The systolic peak is the first hump, not the global maximum.
signal = localTwoHumpPpg(0.85,1.0,0.15,0.32);
result = mpk.PPGProcessor('LabelFrequencyAnalysis',[], ...
    'ReturnIntermediates',true).process(signal);
peakTimes = signal.Time(result.Intermediates.pulse_peak_indices);
verifyEqual(testCase,mod(peakTimes,0.8),repmat(0.15,size(peakTimes)),'AbsTol',0.02)
end

function testCleaningRemovesContentAbove8Hz(testCase)
% A 30 Hz interference must not reach the cleaned signal.
clean = localPpg(repmat(0.8,1,150));
noisy = mpk.Signal(clean.Samples + 0.05*sin(2*pi*30*clean.Time), ...
    clean.Time,100,'ppg','a.u.');
processor = mpk.PPGProcessor('LabelFrequencyAnalysis',[],'ReturnIntermediates',true);
reference = processor.process(clean);
result = processor.process(noisy);
difference = result.Intermediates.cleaned_signal - reference.Intermediates.cleaned_signal;
verifyLessThan(testCase,max(abs(difference(201:end-200))),0.005)
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

function signal = localTwoHumpPpg(first,second,firstAt,secondAt)
% 150 pulses of two Gaussian humps (width 50 ms), one every 0.8 s.
fs = 100;
tau = (0:79)'/fs;
pulse = first*exp(-0.5*((tau - firstAt)/0.05).^2) + ...
    second*exp(-0.5*((tau - secondAt)/0.05).^2);
samples = [repmat(pulse,150,1); 0];
signal = mpk.Signal(samples,(0:numel(samples)-1)'/fs,fs,'ppg','a.u.');
end

function ibis = localRsaIbis(durationSeconds)
ibis = [];
elapsed = 0;
while elapsed < durationSeconds
    ibis(end+1) = 0.9 + 0.05*sin(2*pi*0.25*elapsed); %#ok<AGROW>
    elapsed = elapsed + ibis(end);
end
end
