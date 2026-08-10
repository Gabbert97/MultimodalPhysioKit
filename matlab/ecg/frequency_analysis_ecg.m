function [mean_PLF_final,std_PLF_final,mean_PHF_final,std_PHF_final, ...
          mean_LF_HF_ratio_final,std_LF_HF_ratio_final,frequency_data] = ...
          frequency_analysis_ecg(ibi_tot,length_ibi_phases,condition_columns, ...
          label_frequency_analysis,saving_flag,phase_selection_flag,output_label,gr,varargin)
%FREQUENCY_ANALYSIS_ECG Compute phase-level frequency-domain HRV features.
%
%   [MEAN_PLF,STD_PLF,MEAN_PHF,STD_PHF,MEAN_LF_HF,STD_LF_HF] =
%   FREQUENCY_ANALYSIS_ECG(IBI_TOT,LENGTH_IBI_PHASES,CONDITION_COLUMNS,
%   LABEL_FREQUENCY_ANALYSIS,SAVING_FLAG,PHASE_SELECTION_FLAG,OUTPUT_LABEL,GR)
%   estimates low-frequency (LF), high-frequency (HF), and LF/HF metrics
%   from inter-beat intervals.
%
%   An optional trailing STFT_CONTEXT struct supplies full-recording IBI
%   values, their original timestamps, the recording time range, and an N-by-2
%   matrix of original phase windows. ecg_signal_processing supplies this
%   context in STFT mode. Omitting it preserves compatibility with direct
%   legacy calls.
%
%   Inputs
%   ------
%   IBI_TOT : Inter-beat intervals in seconds. Values are expected to be
%             positive, finite RR/IBI durations.
%   LENGTH_IBI_PHASES : Number of accepted IBI samples assigned to each phase.
%   CONDITION_COLUMNS : Optional phase indices used for condition-specific plots.
%   LABEL_FREQUENCY_ANALYSIS : Analysis mode.
%             1 = STFT on the uniformly resampled tachogram.
%             2 = phase-level PSD / parametric spectral analysis.
%   SAVING_FLAG : If true, diagnostic figures are saved.
%   PHASE_SELECTION_FLAG : Plot-selection flag for condition-specific displays.
%   OUTPUT_LABEL : Prefix used for saved figure filenames.
%   GR : Logical flag. If true, diagnostic figures are generated.
%   STFT_CONTEXT : Optional internal struct with fields IBI, IBITime,
%             TimeRange, and PhaseWindows.
%
%   Outputs
%   -------
%   MEAN_PLF, STD_PLF : Phase-level normalized LF power and its variability.
%   MEAN_PHF, STD_PHF : Phase-level normalized HF power and its variability.
%   MEAN_LF_HF, STD_LF_HF : Phase-level LF/HF ratio and its variability.
%   FREQUENCY_DATA : Scalar struct containing the exact uniformly sampled
%             tachogram and method-specific spectral arrays calculated by
%             this call. In STFT mode, Matrix and Power are
%             numel(Frequencies)-by-numel(Time). In PSD mode, PSD and
%             Frequencies are one-element-per-phase cell arrays.
%
%   Notes
%   -----
%   - Raw IBI sequences are resampled to a uniform 4 Hz tachogram.
%   - LF and HF bands use standard HRV ranges: LF = 0.04-0.15 Hz and
%     HF = 0.15-0.40 Hz.
%   - STFT mode uses a fixed 60-second Hamming window and rejects shorter
%     input phases before this function is called.
%   - With STFT_CONTEXT, one global STFT is mapped to phase windows on the
%     original timeline. A centered STFT window may use surrounding recording
%     samples as context near a phase boundary.
%   - Without STFT_CONTEXT, phase boundaries are derived from IBI counts for
%     backward compatibility and for phase-level PSD analysis.

%% Spectral-analysis configuration
fs_hrv = 4;                  % Uniform tachogram sampling frequency [Hz]
lf_band = [0.04 0.15];       % Low-frequency HRV band [Hz]
hf_band = [0.15 0.40];       % High-frequency HRV band [Hz]
total_band = [0.04 0.40];    % Normalization band [Hz]

winlen_seconds = 60;          % Fixed minimum needed for LF-HRV interpretation
label_PSD = 1;               % 1 Burg, 2 modified covariance, 3 Yule-Walker, 4 Welch

figure_output_folder = 'figures';
if saving_flag && ~isfolder(figure_output_folder)
    mkdir(figure_output_folder)
end

%% Validate and normalize inputs
ibi_tot = ibi_tot(:);
length_ibi_phases = length_ibi_phases(:)';
stft_context = [];
if ~isempty(varargin)
    stft_context = varargin{1};
end

assert(all(isfinite(ibi_tot)) && all(ibi_tot > 0), ...
    'ibi_tot must contain positive finite IBI values in seconds.')

n_phases = numel(length_ibi_phases);
phase_end_counts = cumsum(length_ibi_phases);

assert(phase_end_counts(end) <= numel(ibi_tot), ...
    'length_ibi_phases exceeds the number of available IBIs.')

%% Convert phase lengths from IBI counts to time boundaries
beat_time = cumsum(ibi_tot);
phase_edges_seconds = [0; beat_time(phase_end_counts)];
phase_windows_seconds = [phase_edges_seconds(1:end-1) ...
    phase_edges_seconds(2:end)];

%% Resample the irregular IBI sequence to a uniform tachogram
ibi_time = [0; beat_time];
ibi_values = [ibi_tot(1); ibi_tot];

uniform_time = (0:1/fs_hrv:beat_time(end))';
ibi_uniform = interp1(ibi_time,ibi_values,uniform_time,'pchip');
ibi_uniform = detrend(ibi_uniform,'constant');

if label_frequency_analysis == 1 && ~isempty(stft_context)
    required_fields = {'IBI','IBITime','TimeRange','PhaseWindows'};
    if ~isstruct(stft_context) || ...
            ~all(isfield(stft_context,required_fields))
        error('MultimodalPhysioKit:ECG:InvalidStftContext', ...
            'The global STFT context is incomplete.')
    end

    global_ibi = stft_context.IBI(:);
    global_ibi_time = stft_context.IBITime(:);
    [global_ibi_time,sort_index] = sort(global_ibi_time);
    global_ibi = global_ibi(sort_index);

    recording_time_range = stft_context.TimeRange(:)';
    phase_windows_seconds = stft_context.PhaseWindows;
    if size(phase_windows_seconds,1) ~= n_phases || ...
            size(phase_windows_seconds,2) ~= 2
        error('MultimodalPhysioKit:ECG:InvalidStftContext', ...
            'Global STFT phase windows must be N-by-2.')
    end

    interpolation_time = [recording_time_range(1); global_ibi_time; ...
        recording_time_range(2)];
    interpolation_values = [global_ibi(1); global_ibi; global_ibi(end)];
    [interpolation_time,unique_index] = unique(interpolation_time,'stable');
    interpolation_values = interpolation_values(unique_index);

    uniform_time = (recording_time_range(1):1/fs_hrv: ...
        recording_time_range(2))';
    ibi_uniform = interp1(interpolation_time,interpolation_values, ...
        uniform_time,'pchip');
    ibi_uniform = detrend(ibi_uniform,'constant');
end

if label_frequency_analysis == 1
    winlen = round(winlen_seconds*fs_hrv);
    if numel(ibi_uniform) < winlen
        error('MultimodalPhysioKit:ECG:StftRecordingTooShort', ...
            'The resampled tachogram is shorter than the fixed 60-second STFT window.')
    end
else
    winlen = max(16,min(round(winlen_seconds*fs_hrv),numel(ibi_uniform)));
end
nfft = max(256,2^nextpow2(winlen));

mean_PLF = nan(1,n_phases);
std_PLF = nan(1,n_phases);
mean_PHF = nan(1,n_phases);
std_PHF = nan(1,n_phases);
mean_LF_HF_ratio = nan(1,n_phases);
std_LF_HF_ratio = nan(1,n_phases);
psd_by_phase = cell(1,n_phases);
psd_frequencies_by_phase = cell(1,n_phases);
psd_lf_power = nan(1,n_phases);
psd_hf_power = nan(1,n_phases);
psd_total_power = nan(1,n_phases);

%% Select original phase markers for optional plots
plot_phase_windows = phase_windows_seconds;

if phase_selection_flag ~= 0 && ~isempty(condition_columns)
    selected_phases = condition_columns(:)';
    selected_phases = selected_phases(selected_phases >= 1 & selected_phases <= n_phases);

    if ~isempty(selected_phases)
        plot_phase_windows = phase_windows_seconds(selected_phases,:);
    end
end

plot_phase_starts = unique(plot_phase_windows(:,1));
plot_phase_ends = unique(plot_phase_windows(:,2));
plot_phase_edges = unique(plot_phase_windows(:));
plot_time_limits = [min(plot_phase_starts) max(plot_phase_ends)];

%% STFT-based HRV analysis
if label_frequency_analysis == 1

    fprintf('The NFFT for the STFT is: %d\n',nfft)

    [S,t_stft,f] = STFT(ibi_uniform,nfft,fs_hrv,winlen,hamming(winlen));
    t_stft = t_stft + uniform_time(1);
    P = abs(S).^2;

    lf_index = f >= lf_band(1) & f <= lf_band(2);
    hf_index = f >  hf_band(1) & f <= hf_band(2);
    total_index = f >= total_band(1) & f <= total_band(2);

    lf_power = sum(P(lf_index,:),1);
    hf_power = sum(P(hf_index,:),1);
    total_power = sum(P(total_index,:),1);

    PLF = lf_power ./ max(total_power,eps);
    PHF = hf_power ./ max(total_power,eps);
    LF_HF_ratio = lf_power ./ max(hf_power,eps);

    for i_phase = 1:n_phases
        phase_index = t_stft >= phase_windows_seconds(i_phase,1) & ...
                      t_stft <= phase_windows_seconds(i_phase,2);

        mean_PLF(i_phase) = mean(PLF(phase_index),'omitnan');
        std_PLF(i_phase) = std(PLF(phase_index),'omitnan');

        mean_PHF(i_phase) = mean(PHF(phase_index),'omitnan');
        std_PHF(i_phase) = std(PHF(phase_index),'omitnan');

        mean_LF_HF_ratio(i_phase) = mean(LF_HF_ratio(phase_index),'omitnan');
        std_LF_HF_ratio(i_phase) = std(LF_HF_ratio(phase_index),'omitnan');
    end

    if gr
        freq_plot_index = f >= 0.03 & f <= 0.50;

        figure
        mesh(t_stft,f(freq_plot_index),P(freq_plot_index,:),'FaceColor','interp')
        view(2)
        hold on
        % Draw solid red phase ends first, then dashed white phase starts so
        % coincident boundaries remain distinguishable by color and style.
        end_lines = xline(plot_phase_ends,'-','Color',[1 0 0], ...
            'LineWidth',1.2);
        start_lines = xline(plot_phase_starts,'--','Color',[1 1 1], ...
            'LineWidth',1.2);
        set(end_lines,'HandleVisibility','off')
        set(start_lines,'HandleVisibility','off')
        set(end_lines(1),'HandleVisibility','on','DisplayName','Phase end')
        set(start_lines(1),'HandleVisibility','on','DisplayName','Phase start')
        yline(0.15,'w','LF/HF boundary')
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        title(sprintf('STFT HRV power, window = %.1f s',winlen/fs_hrv))
        ylim([0.04 0.40])
        xlim(plot_time_limits)
        grid on
        legend([start_lines(1) end_lines(1)], ...
            {'Phase start','Phase end'},'TextColor','white', ...
            'Color','none','Location','best')
        set(gca,'FontSize',15)

        if saving_flag
            saveas(gcf,fullfile(figure_output_folder, ...
                strcat(output_label,'_STFT_HRV.png')))
        end

        figure
        plot(t_stft,LF_HF_ratio,'LineWidth',1.5)
        hold on
        xline(plot_phase_edges,'k--','LineWidth',0.8)
        xlabel('Time (s)')
        ylabel('LF/HF (-)')
        title('STFT LF/HF')
        xlim(plot_time_limits)
        grid on
        set(gca,'FontSize',15)

        if saving_flag
            saveas(gcf,fullfile(figure_output_folder, ...
                strcat(output_label,'_STFT_LFHF.png')))
        end
    end

%% Phase-level PSD / parametric HRV analysis
elseif label_frequency_analysis == 2

    for i_phase = 1:n_phases
        phase_index = uniform_time >= phase_edges_seconds(i_phase) & ...
                      uniform_time <= phase_edges_seconds(i_phase + 1);

        x_phase = ibi_uniform(phase_index);

        if numel(x_phase) < 16
            warning('Phase %d is too short for HRV frequency analysis.',i_phase)
            continue
        end

        [pxx,f] = local_psd(x_phase,fs_hrv,nfft,label_PSD);

        lf_power = local_band_power(pxx,f,lf_band);
        hf_power = local_band_power(pxx,f,hf_band);
        total_power = local_band_power(pxx,f,total_band);
        psd_by_phase{i_phase} = pxx;
        psd_frequencies_by_phase{i_phase} = f;
        psd_lf_power(i_phase) = lf_power;
        psd_hf_power(i_phase) = hf_power;
        psd_total_power(i_phase) = total_power;

        mean_PLF(i_phase) = lf_power / max(total_power,eps);
        mean_PHF(i_phase) = hf_power / max(total_power,eps);
        mean_LF_HF_ratio(i_phase) = lf_power / max(hf_power,eps);

        std_PLF(i_phase) = 0;
        std_PHF(i_phase) = 0;
        std_LF_HF_ratio(i_phase) = 0;
    end

    if gr
        figure
        plot(1:n_phases,mean_LF_HF_ratio,'o-','LineWidth',1.5)
        xlabel('Phase (-)')
        ylabel('LF/HF (-)')
        title('Phase-level LF/HF')
        grid on
        set(gca,'FontSize',15)

        if saving_flag
            saveas(gcf,fullfile(figure_output_folder, ...
                strcat(output_label,'_PSD_LFHF.png')))
        end
    end

else
    error(['Invalid label_frequency_analysis. Use:\n' ...
           '1 = STFT analysis\n' ...
           '2 = phase-level PSD analysis'])
end

%% Return phase-level features
mean_PLF_final = mean_PLF;
std_PLF_final = std_PLF;
mean_PHF_final = mean_PHF;
std_PHF_final = std_PHF;
mean_LF_HF_ratio_final = mean_LF_HF_ratio;
std_LF_HF_ratio_final = std_LF_HF_ratio;

frequency_data = struct( ...
    'Method',label_frequency_analysis, ...
    'SamplingFrequency',fs_hrv, ...
    'UniformTime',uniform_time, ...
    'UniformIBI',ibi_uniform, ...
    'PhaseWindows',phase_windows_seconds, ...
    'LFBand',lf_band, ...
    'HFBand',hf_band, ...
    'TotalBand',total_band, ...
    'NFFT',nfft);

if label_frequency_analysis == 1
    frequency_data.Matrix = S;
    frequency_data.Power = P;
    frequency_data.Frequencies = f;
    frequency_data.Time = t_stft;
    frequency_data.LFPower = lf_power;
    frequency_data.HFPower = hf_power;
    frequency_data.TotalPower = total_power;
    frequency_data.PLF = PLF;
    frequency_data.PHF = PHF;
    frequency_data.LFHF = LF_HF_ratio;
    frequency_data.WindowLengthSamples = winlen;
    frequency_data.WindowLengthSeconds = winlen_seconds;
else
    frequency_data.PSD = psd_by_phase;
    frequency_data.Frequencies = psd_frequencies_by_phase;
    frequency_data.LFPower = psd_lf_power;
    frequency_data.HFPower = psd_hf_power;
    frequency_data.TotalPower = psd_total_power;
    frequency_data.PLF = mean_PLF;
    frequency_data.PHF = mean_PHF;
    frequency_data.LFHF = mean_LF_HF_ratio;
end

end

%% Local spectral helper functions
function [pxx,f] = local_psd(x,fs,nfft,label_PSD)
%LOCAL_PSD Estimate one-sided PSD using the requested method.
% Falls back to a basic FFT-based PSD if toolbox functions are unavailable.

x = detrend(x(:),'constant');

switch label_PSD
    case 1
        order = max(4,min(16,floor(numel(x)/4)));
        if exist('pburg','file')
            [pxx,f] = pburg(x,order,nfft,fs);
        else
            [pxx,f] = local_fft_psd(x,fs,nfft);
        end

    case 2
        order = max(4,min(16,floor(numel(x)/4)));
        if exist('pmcov','file')
            [pxx,f] = pmcov(x,order,nfft,fs);
        else
            [pxx,f] = local_fft_psd(x,fs,nfft);
        end

    case 3
        order = max(4,min(16,floor(numel(x)/4)));
        if exist('pyulear','file')
            [pxx,f] = pyulear(x,order,nfft,fs);
        else
            [pxx,f] = local_fft_psd(x,fs,nfft);
        end

    case 4
        if exist('pwelch','file')
            [pxx,f] = pwelch(x,[],[],nfft,fs);
        else
            [pxx,f] = local_fft_psd(x,fs,nfft);
        end

    otherwise
        [pxx,f] = local_fft_psd(x,fs,nfft);
end

end

function power_value = local_band_power(pxx,f,band)
%LOCAL_BAND_POWER Integrate PSD values within a frequency band.

index = f >= band(1) & f <= band(2);

if ~any(index)
    power_value = NaN;
    return
end

df = mean(diff(f));
power_value = sum(pxx(index)) * df;

end

function [pxx,f] = local_fft_psd(x,fs,nfft)
%LOCAL_FFT_PSD Basic one-sided FFT PSD fallback.

x = x(:);
w = hamming(numel(x));
xw = x .* w;

X = fft(xw,nfft);
n_freq = floor(nfft/2) + 1;

pxx = abs(X(1:n_freq)).^2 / (fs * sum(w.^2));
f = (0:n_freq-1)' * fs / nfft;

end
