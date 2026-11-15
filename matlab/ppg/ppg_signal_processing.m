function [features,intermediates] = ppg_signal_processing( ...
        ppg,time,fs,label_frequency_analysis,show_plots)
%PPG_SIGNAL_PROCESSING Detect PPG pulses and extract features from one signal.
%
% Inputs
% ------
% PPG   : PPG samples (any amplitude unit).
% TIME  : Strictly increasing time vector [s], same length as PPG.
% FS    : Sampling frequency [Hz]; must exceed 6 Hz.
% LABEL_FREQUENCY_ANALYSIS : 1 = STFT, 2 = PSD (both need >= 60 s), or
%         [] to skip HRV frequency analysis (e.g. short real-time windows).
% SHOW_PLOTS : Logical; if true, diagnostic figures are drawn.
%
% Outputs
% -------
% FEATURES      : N-by-2 cell array {feature_name, value}.
% INTERMEDIATES : Struct with the cleaned and filtered signals, pulse feet/peaks
%                 (one-based indexes), per-pulse values and validity, and
%                 accepted IBIs.
%
% The supplied signal is analysed as one piece (recording, phase, or a
% completed real-time window); all thresholds adapt to the pulses found.
% Beats are located on the 0.5-3.0 Hz band-passed signal; all pulse
% characteristics are measured on the cleaned signal (0.5 Hz high-pass only),
% which removes baseline drift but keeps the pulse morphology.
% Pulse k spans feet k..k+1. A pulse is valid when its duration is within
% 0.8-1.2x the running median of up to 11 neighbouring pulses and its
% systolic peak lies in the first half. IBIs (peak to peak) are kept only
% between two consecutive valid pulses with 40-180 bpm; rejected intervals
% are skipped, never imputed.

if nargin < 5
    show_plots = false;
end
ppg = double(ppg(:));
time = double(time(:));
maxBpm = 180;
minBpm = 40;

if 3.0 >= fs/2
    error('MultimodalPhysioKit:PPG:InvalidSamplingFrequency', ...
        'Sampling frequency must exceed 6 Hz for 0.5-3.0 Hz PPG filtering.')
end
if ~isempty(label_frequency_analysis) && numel(ppg)/fs < 60
    error('MultimodalPhysioKit:PPG:TooShortForFrequencyAnalysis', ...
        ['PPG frequency analysis requires at least 60 seconds; use an ' ...
         'empty LabelFrequencyAnalysis for shorter signals.'])
end

%% Pulse detection
cleaned = filtering(ppg,fs,0.5,0,3,0,0);             % drift removed, shape kept
filtered = filtering(cleaned,fs,0,3.0,0,3,0);        % 0.5-3.0 Hz: 30-180 bpm
[~,troughs] = findpeaks(-filtered, ...
    'MinPeakDistance',max(1,round(fs*60/maxBpm)));
if numel(troughs) < 3
    error('MultimodalPhysioKit:PPG:TooFewPulses', ...
        'Too few PPG pulses were detected.')
end
% Move each filtered trough to the cleaned minimum within -10%/+30% of the
% median trough distance: these are the pulse feet (onsets).
medianDistance = median(diff(troughs));
before = round(0.1*medianDistance);
after = round(0.3*medianDistance);
feet = zeros(numel(troughs),1);
for k = 1:numel(troughs)
    lo = max(1,troughs(k) - before);
    hi = min(numel(cleaned),troughs(k) + after);
    [~,offset] = min(cleaned(lo:hi));
    feet(k) = lo + offset - 1;
end
feet = unique(feet);
if numel(feet) < 3
    error('MultimodalPhysioKit:PPG:TooFewPulses', ...
        'Too few PPG pulses were detected.')
end

starts = feet(1:end-1);
stops = feet(2:end);
peaks = zeros(size(starts));
for k = 1:numel(starts)
    [~,offset] = max(cleaned(starts(k):stops(k)));
    peaks(k) = starts(k) + offset - 1;
end
durations = time(stops) - time(starts);
riseTimes = time(peaks) - time(starts);
baselineAtPeak = cleaned(starts) + (cleaned(stops) - cleaned(starts)) .* riseTimes ./ durations;
amplitudes = cleaned(peaks) - baselineAtPeak;

reference = zeros(size(durations));
for k = 1:numel(durations)
    reference(k) = median(durations(max(1,k-5):min(numel(durations),k+5)));
end
valid = durations >= 0.8*reference & durations <= 1.2*reference & ...
    peaks > starts & peaks < stops & (peaks - starts) <= (stops - starts)/2;

%% Accepted inter-beat intervals
peakTimes = time(peaks);
intervals = diff(peakTimes);
bpmAll = 60 ./ intervals;
flags = valid(1:end-1) & valid(2:end) & bpmAll >= minBpm & bpmAll <= maxBpm;
ibi = intervals(flags);
ibiStartTimes = peakTimes(find(flags));        %#ok<FNDSB>
ibiEndTimes = peakTimes(find(flags) + 1);
if nnz(valid) < 3 || numel(ibi) < 3
    error('MultimodalPhysioKit:PPG:TooFewPulses', ...
        'Too few valid PPG pulses were detected.')
end
bpm = 60 ./ ibi;

%% Features
names = {};
values = {};
pulseValues = {amplitudes(valid),durations(valid),riseTimes(valid)};
pulseNames = {'ppg_pulse_amplitude','ppg_pulse_duration','ppg_rise_time'};
for k = 1:numel(pulseNames)
    names = [names, strcat(pulseNames{k},{'_mean','_std','_median'})]; %#ok<AGROW>
    values = [values, {mean(pulseValues{k}),std(pulseValues{k}), ...
        median(pulseValues{k})}]; %#ok<AGROW>
end
[pnn50,pnn50Max] = ecg_pnn50_features(ibi,flags);
names = [names, {'ppg_bpm_mean','ppg_bpm_std','ppg_bpm_median', ...
    'ppg_pnn50','ppg_pnn50_max'}];
values = [values, {mean(bpm),std(bpm),median(bpm),pnn50,pnn50Max}];

frequencyData = struct();
if ~isempty(label_frequency_analysis)
    stftContext = [];
    if label_frequency_analysis == 1
        stftContext = struct( ...
            'IBI',ibi, ...
            'IBITime',(ibiStartTimes + ibiEndTimes)/2, ...
            'TimeRange',[time(1) time(end)], ...
            'PhaseWindows',[time(1) time(end)]);
    end
    [plf,plfStd,phf,phfStd,lfHf,lfHfStd,frequencyData] = ...
        frequency_analysis_ecg(ibi,numel(ibi),1,label_frequency_analysis, ...
            false,0,'ppg',show_plots,stftContext,{ibiStartTimes});
    names = [names, {'ppg_plf_mean','ppg_plf_std','ppg_phf_mean', ...
        'ppg_phf_std','ppg_lf_hf_mean','ppg_lf_hf_std'}];
    values = [values, {plf,plfStd,phf,phfStd,lfHf,lfHfStd}];
end
features = [names(:), values(:)];

intermediates = struct( ...
    'cleaned_signal',cleaned(:), ...
    'filtered_signal',filtered(:), ...
    'pulse_onset_indices',feet(:), ...
    'pulse_peak_indices',peaks(:), ...
    'pulse_peak_times',peakTimes(:), ...
    'pulse_valid',valid(:), ...
    'pulse_amplitudes',amplitudes(:), ...
    'pulse_durations',durations(:), ...
    'pulse_rise_times',riseTimes(:), ...
    'ibi',ibi(:), ...
    'ibi_times',ibiStartTimes(:), ...
    'bpm',bpm(:), ...
    'frequency',frequencyData);

%% Diagnostic plots
if show_plots
    figure('Name','PPG pulses')
    subplot(3,1,1)
    plot(time,ppg,'Color',[0.4 0.4 0.4])
    ylabel('Raw PPG')
    title(sprintf('PPG pulses: %d accepted, %d rejected', ...
        nnz(valid),nnz(~valid)))
    grid on
    subplot(3,1,2)
    hold on
    for k = 1:numel(starts)
        segment = starts(k):stops(k);
        if valid(k)
            plot(time(segment),cleaned(segment),'b')
        else
            plot(time(segment),cleaned(segment),'r')
        end
    end
    plot(time(feet),cleaned(feet),'kv')
    plot(time(peaks(valid)),cleaned(peaks(valid)),'bo')
    plot(time(peaks(~valid)),cleaned(peaks(~valid)),'rx')
    ylabel('Cleaned PPG (>0.5 Hz)')
    grid on
    subplot(3,1,3)
    plot(time,filtered)
    xlabel('Time (s)')
    ylabel('Band-pass 0.5-3.0 Hz')
    grid on

    figure('Name','PPG BPM trend')
    plot(ibiStartTimes,bpm,'o-','LineWidth',1.2)
    xlabel('Time (s)')
    ylabel('BPM')
    title('PPG BPM trend (accepted IBIs)')
    grid on
end
end
