function [ecg_features,intermediates] = ecg_signal_processing(ecg,time,i_phases,condition_columns,fs,label_frequency_analysis,saving_flag,phase_selection_flag,output_label,gr)
%ECG_SIGNAL_PROCESSING Extract phase-level ECG, HR, and HRV features.
%
%   ECG_FEATURES = ECG_SIGNAL_PROCESSING(ECG,TIME,I_PHASES,CONDITION_COLUMNS,FS,
%   LABEL_FREQUENCY_ANALYSIS,SAVING_FLAG,PHASE_SELECTION_FLAG,OUTPUT_LABEL,GR)
%   processes a converted ECG signal, detects R-peaks independently in each
%   phase for time-domain features, and appends frequency-domain HRV metrics.
%   STFT mode additionally detects R-peaks once across the complete recording
%   to construct a nonduplicated tachogram on the original time axis.
%
%   Inputs
%   ------
%   ECG    : Converted ECG vector.
%   TIME   : Time vector in seconds, same length as ECG.
%   I_PHASES : 2-by-N matrix of phase limits in samples. Row 1 contains
%              phase starts; row 2 contains phase ends; columns are phases.
%   CONDITION_COLUMNS : Phase indices used for condition-specific frequency
%              plots when PHASE_SELECTION_FLAG requests them.
%   FS     : ECG sampling frequency, in Hz.
%   LABEL_FREQUENCY_ANALYSIS : Frequency-analysis mode.
%              1 = STFT-based HRV analysis.
%              2 = phase-level PSD / parametric HRV analysis.
%   SAVING_FLAG : If true, generated figures are saved to disk.
%   PHASE_SELECTION_FLAG : Plot-selection flag forwarded to frequency analysis.
%   OUTPUT_LABEL : Prefix used when saving diagnostic figures.
%   GR     : Logical flag. If true, diagnostic figures are generated.
%
%   Output
%   ------
%   ECG_FEATURES : Cell array whose first column contains feature names and
%                  whose remaining columns contain phase-level feature values.
%   INTERMEDIATES : Scalar struct containing phase-local detector, screened
%                  IBI/BPM, and frequency-analysis values already calculated
%                  by the scientific pipeline.
%
%   Notes
%   -----
%   - Time-domain R-peak screening state is reset for every phase. Overlapping
%     samples are therefore reevaluated independently in overlapping phases.
%   - STFT mode uses one full-recording tachogram and maps spectral values to
%     the original phase windows. It requires phases of at least 60 seconds.
%   - Phase-level PSD mode retains the accepted phase-specific IBI workflow.
%   - A one-phase call on an extracted Signal calculates its own complete
%     tachogram and its own STFT or PSD from that supplied Signal only.
%   - The function expects helper functions from the project MATLAB folder.

fprintf('ECG HR and HRV processing in progress. Please wait..\n')

if label_frequency_analysis == 1
    phase_durations_seconds = ...
        (i_phases(2,:) - i_phases(1,:) + 1) / fs;
    if any(phase_durations_seconds < 60)
        error('MultimodalPhysioKit:ECG:StftWindowTooShort', ...
            ['STFT HRV analysis requires every input window to be at ' ...
             'least 60 seconds. Shorter windows do not support ' ...
             'meaningful estimation of the LF band (0.04-0.15 Hz).'])
    end
end

%% Configuration and initialization
check_flag = gr;
up_percent_threshold_amp = 99;
down_percent_treshold_amp = 1;
default_l_window_amp = 10;                                                      % [-] number of peaks
l_window_time = 5;

indexes_phases = split_signal_phases(i_phases);
number_phases = length(indexes_phases);

BPM_def = [];
length_BPM_phases = [];
ibi_def = [];

ecg_features = {'ecg_bpm_mean';'ecg_bpm_std';'ecg_bpm_median';'ecg_pnn50';'ecg_pnn50_max';'ecg_plf_mean';'ecg_plf_std';'ecg_phf_mean';'ecg_phf_std';'ecg_lf_hf_mean';'ecg_lf_hf_std'};
intermediates = struct( ...
    'filtered_signal',{cell(1,number_phases)}, ...
    'candidate_r_peak_indices',{cell(1,number_phases)}, ...
    'candidate_r_peak_amplitudes',{cell(1,number_phases)}, ...
    'accepted_r_peak_indices',{cell(1,number_phases)}, ...
    'accepted_r_peak_amplitudes',{cell(1,number_phases)}, ...
    'ibi',{cell(1,number_phases)}, ...
    'ibi_times',{cell(1,number_phases)}, ...
    'bpm',{cell(1,number_phases)}, ...
    'bpm_times',{cell(1,number_phases)}, ...
    'frequency',struct());

for phase = 1:number_phases

    % Keep adaptive peak-screening state independent between phases.
    l_window_amp = default_l_window_amp;
    ecg_phase = ecg(indexes_phases{phase});
    time_phase = time(indexes_phases{phase});

    %% R-peak detection
    % Detect candidate R-peaks on the filtered ECG returned by
    % My_R_peaks_detection.
    [qrs_amp_raw_phase,qrs_i_raw_phase,filtered_ecg_phase] = My_R_peaks_detection(ecg_phase,time_phase,fs,0);
    intermediates.filtered_signal{phase} = filtered_ecg_phase;
    intermediates.candidate_r_peak_indices{phase} = qrs_i_raw_phase;
    intermediates.candidate_r_peak_amplitudes{phase} = ...
        qrs_amp_raw_phase;
    
    %% R-peak amplitude screening
    % Define upper and lower amplitude thresholds from the detected
    % R-peak amplitude distribution.
    %
    up_initial_threshold_amp_r_peak = prctile(qrs_amp_raw_phase,up_percent_threshold_amp);
    down_initial_threshold_amp_r_peak = prctile(qrs_amp_raw_phase,down_percent_treshold_amp);
 

    % Use the initial accepted peaks to seed the adaptive amplitude check.
    % Candidate amplitudes are compared against percentile-based limits:
    %
    % 1) upper threshold: UP_PERCENT_THRESHOLD_AMP percentile.
    % 2) lower threshold: DOWN_PERCENT_THRESHOLD_AMP percentile.
    qrs_i_initial = qrs_i_raw_phase(qrs_amp_raw_phase < up_initial_threshold_amp_r_peak & qrs_amp_raw_phase > down_initial_threshold_amp_r_peak);
    qrs_amp_initial = qrs_amp_raw_phase(qrs_amp_raw_phase < up_initial_threshold_amp_r_peak & qrs_amp_raw_phase > down_initial_threshold_amp_r_peak);
    
    if length(qrs_i_initial) < l_window_amp
        l_window_amp = round(l_window_amp/2);
    end
    
    qrs_i_def = qrs_i_initial(1:l_window_amp);
    qrs_amp_def = qrs_amp_initial(1:l_window_amp);
    count = find(qrs_i_raw_phase == qrs_i_def(end));

    % Continue accepting peaks whose amplitude remains consistent with the
    n_qrs_peaks = length(qrs_i_raw_phase) - count;
    for i = 1 : n_qrs_peaks
        mean_amp_qrs_peaks = mean(qrs_amp_def(length(qrs_amp_def)-l_window_amp+1:end));
        % recent accepted R-peak amplitude range.
        %
        if qrs_amp_raw_phase(count + i) >= mean_amp_qrs_peaks * 0.6 && qrs_amp_raw_phase(count + i) <= mean_amp_qrs_peaks * 1.5
            qrs_i_def = [qrs_i_def qrs_i_raw_phase(count + i)];
            qrs_amp_def = [qrs_amp_def qrs_amp_raw_phase(count + i)];
        end
    end
    intermediates.accepted_r_peak_indices{phase} = qrs_i_def;
    intermediates.accepted_r_peak_amplitudes{phase} = qrs_amp_def;

    if check_flag
        figure
        fig_rr_plot(1) = subplot(2,1,1);
        plot(time_phase,filtered_ecg_phase)
        hold on
        for i = 1:length(qrs_i_raw_phase)
            plot(time_phase(qrs_i_raw_phase),filtered_ecg_phase(qrs_i_raw_phase),'*r')
        end
        for i = 1:length(qrs_i_def)
            plot(time_phase(qrs_i_def),filtered_ecg_phase(qrs_i_def),'*b')
        end
        xlabel('time (s)')
        ylabel('ECG (mV)')
        grid on
        xlim([time_phase(1) time_phase(end)])
        set(gca,'FontSize',14)
        fig_rr_plot(2) = subplot(2,1,2);
        hold on
        for i = 1:length(qrs_i_raw_phase)
            xx = ones(2,1)*time_phase(qrs_i_raw_phase(i));
            yy = [0 qrs_amp_raw_phase(i)];
            plot(xx,yy,'-r')
        end
        for i = 1:length(qrs_i_def)
            xx = ones(2,1)*time_phase(qrs_i_def(i));
            yy = [0 qrs_amp_def(i)];
            plot(xx,yy,'-b')
        end
        xlabel('time (s)')
        ylabel('Amplitude ECG (mV)')
        grid on
        xlim([time_phase(1) time_phase(end)])
        linkaxes(fig_rr_plot,'x')
        set(gca,'FontSize',14)
    end

    %% IBI validation
    % Convert accepted R-peak intervals from samples to seconds.
    diff_i_r_peaks = diff(qrs_i_def);
    diff_t_r_peaks = diff_i_r_peaks/fs;

    control_ibi_flag = zeros(length(diff_t_r_peaks),1);

    mean_diff_t_peaks = mean(diff_t_r_peaks);

    ibi_preliminary = diff_t_r_peaks(diff_t_r_peaks < 1.5 * mean_diff_t_peaks & diff_t_r_peaks > 0.5 * mean_diff_t_peaks);
    indexes_ibi_preliminary = qrs_i_def(diff_t_r_peaks < 1.5 * mean_diff_t_peaks & diff_t_r_peaks > 0.5 * mean_diff_t_peaks);

    ibi_def_phase = ibi_preliminary(1:l_window_time);
    indexes_ibi_def = indexes_ibi_preliminary(1:l_window_time);

    for i = 1:l_window_time
        control_ibi_flag(find(qrs_i_def == indexes_ibi_def(i))) = 1;
    end

    count_ibi = find(qrs_i_def == indexes_ibi_def(end));
    n_ibi = length(diff_t_r_peaks) - count_ibi;
    for i = 1 : n_ibi
        mean_diff_t_peaks =  mean(ibi_def_phase(length(ibi_def_phase)-l_window_time+1:end));
        if (diff_t_r_peaks(count_ibi + i) >= mean_diff_t_peaks * 0.7) && (diff_t_r_peaks(count_ibi + i) <= mean_diff_t_peaks * 1.5)
            ibi_def_phase = [ibi_def_phase diff_t_r_peaks(count_ibi + i)];
            control_ibi_flag(count_ibi + i) = 1;
        end
    end

    if check_flag
        j = 1;
        fig_rr_plot(1) = subplot(2,1,1);
        for i = 1 : length(qrs_i_def)-1
            if control_ibi_flag(i) == 1 && qrs_i_def(i) ~= qrs_i_def(end)
                plot([time_phase(qrs_i_def(i)) (time_phase(qrs_i_def(i))+ibi_def_phase(j))],[0 0],'-k','LineWidth',1.5)
                plot(time_phase(qrs_i_def(i))+ibi_def_phase(j),0,'*k')
                j = j + 1;
            end
        end
    end

    %% Time-domain ECG and HRV features
    % Convert inter-beat intervals to heart rate in beats per minute.
    raw_BPM = ibi_def_phase;
    raw_BPM = 60./raw_BPM;
    % Replace physiologically implausible BPM values using neighboring samples.
    i_BPM_outlier = find(raw_BPM < 40);
    i_BPM_outlier = [i_BPM_outlier;find(raw_BPM > 160)];
    BPM = raw_BPM;
    
    % Case 1: outliers are internal BPM samples.
    if ~isempty(i_BPM_outlier) && ~any(i_BPM_outlier == 1) && ~any(i_BPM_outlier == length(raw_BPM))
        BPM(i_BPM_outlier) = (raw_BPM(i_BPM_outlier + 1) + raw_BPM(i_BPM_outlier - 1))/2;
    % Case 2: outliers include both the first and last BPM samples.
    elseif ~isempty(i_BPM_outlier) && any(i_BPM_outlier == 1) && any(i_BPM_outlier == length(raw_BPM)) 
        BPM(i_BPM_outlier(1)) = (raw_BPM(i_BPM_outlier(1) + 1) + raw_BPM(i_BPM_outlier(1) + 2))/2;
        BPM(i_BPM_outlier(end)) = (raw_BPM(i_BPM_outlier(end) - 1) + raw_BPM(i_BPM_outlier(end) - 2))/2;
        len_ok_i_outlier = length(i_BPM_outlier(2:end-1));
        if ~isempty(len_ok_i_outlier)
            BPM(i_BPM_outlier(2:end-1)) = (raw_BPM(i_BPM_outlier(2:end-1) - 1) + raw_BPM(i_BPM_outlier(2:end-1) + 1))/2;
        end
    elseif  ~isempty(i_BPM_outlier) && any(i_BPM_outlier == 1)
        BPM(i_BPM_outlier(1)) = (raw_BPM(i_BPM_outlier(1) + 1) + raw_BPM(i_BPM_outlier(1) + 2))/2;
        len_ok_i_outlier = length(i_BPM_outlier(2:end));
        if ~isempty(len_ok_i_outlier)
            BPM(i_BPM_outlier(2:end)) = (raw_BPM(i_BPM_outlier(2:end) - 1) + raw_BPM(i_BPM_outlier(2:end) + 1))/2;
        end
    elseif  ~isempty(i_BPM_outlier) && any(i_BPM_outlier == length(raw_BPM))
        BPM(i_BPM_outlier(end)) = (raw_BPM(i_BPM_outlier(end) - 1) + raw_BPM(i_BPM_outlier(end) - 2))/2;
        len_ok_i_outlier = length(i_BPM_outlier(2:end));
        if ~isempty(len_ok_i_outlier)
            BPM(i_BPM_outlier(1:end-1)) = (raw_BPM(i_BPM_outlier(1:end-1) - 1) + raw_BPM(i_BPM_outlier(1:end-1) + 1))/2;
        end
    end

    accepted_ibi_positions = find(control_ibi_flag ~= 0);
    accepted_ibi_times = time_phase(qrs_i_def(accepted_ibi_positions));
    intermediates.ibi{phase} = ibi_def_phase;
    intermediates.ibi_times{phase} = accepted_ibi_times;
    intermediates.bpm{phase} = BPM;
    intermediates.bpm_times{phase} = accepted_ibi_times;

    %% Time-domain HRV features
    % Mean, standard deviation, and median heart rate.
    ecg_features{1,phase+1} = mean(BPM);
    ecg_features{2,phase+1} = std(BPM);
    ecg_features{3,phase+1} = median(BPM);
    % pNN50 (share of successive differences > 50 ms between adjacent
    % accepted IBIs) and the longest run of such differences.
    [ecg_features{4,phase+1},ecg_features{5,phase+1}] = ...
        ecg_pnn50_features(ibi_def_phase,control_ibi_flag);
    length_BPM_phases = [length_BPM_phases; length(BPM)];
    BPM_def = [BPM_def,BPM];
    ibi_def = [ibi_def,ibi_def_phase];
end

if gr
    x_line_indexes = cumsum(length_BPM_phases);
    figure
    hold on
    phase_start_index = 1;
    phase_colors = lines(number_phases);
    for phase = 1:number_phases
        phase_end_index = x_line_indexes(phase);
        phase_plot_indexes = phase_start_index:phase_end_index;
        plot(phase_plot_indexes,BPM_def(phase_plot_indexes), ...
            'Color',phase_colors(phase,:),'LineWidth',1.5)
        phase_start_index = phase_end_index + 1;
    end
    xline(x_line_indexes,'--')
    grid on
    xlabel('Accepted IBI sequence by phase (overlap may repeat beats)')
    ylabel('BPM')
    title('Phase-wise BPM trend')

    y_limits = ylim;
    label_y = y_limits(2) - 0.03*diff(y_limits);
    phase_start_indexes = [1; x_line_indexes(1:end-1) + 1];
    for phase = 1:number_phases
        label_x = (phase_start_indexes(phase) + ...
            x_line_indexes(phase)) / 2;
        text(label_x,label_y,sprintf('Phase %d',phase), ...
            'HorizontalAlignment','center','VerticalAlignment','top', ...
            'Rotation',90,'FontSize',9)
    end
    set(gca,'FontSize',15)
end

%% Frequency-domain HRV features
stft_context = [];
if label_frequency_analysis == 1
    stft_context = local_build_stft_context(ecg,time,i_phases,fs);
end

[mean_PLF_final,std_PLF_final,mean_PHF_final,std_PHF_final,mean_LF_HF_ratio_final,std_LF_HF_ratio_final,frequency_data] = frequency_analysis_ecg(ibi_def,length_BPM_phases,condition_columns,label_frequency_analysis,saving_flag,phase_selection_flag,output_label,gr,stft_context,intermediates.ibi_times);
intermediates.frequency = frequency_data;
for phase = 1 : number_phases
    ecg_features{6,phase+1} = mean_PLF_final(phase);
    ecg_features{7,phase+1} = std_PLF_final(phase);
    ecg_features{8,phase+1} = mean_PHF_final(phase);
    ecg_features{9,phase+1} = std_PHF_final(phase);
    ecg_features{10,phase+1} = mean_LF_HF_ratio_final(phase);
    ecg_features{11,phase+1} = std_LF_HF_ratio_final(phase);
end
end

function context = local_build_stft_context(ecg,time,i_phases,fs)
%LOCAL_BUILD_STFT_CONTEXT Build a full-recording tachogram for STFT analysis.
% R-peaks are detected once on the original ECG. IBI timestamps and phase
% windows retain the recording time axis, including overlapping windows.

[~,qrs_indexes] = My_R_peaks_detection(ecg,time,fs,0);
qrs_indexes = qrs_indexes(:);
if numel(qrs_indexes) < 3
    error('MultimodalPhysioKit:ECG:InsufficientGlobalRPeaks', ...
        'At least three R peaks are required for global STFT analysis.')
end

ibi = diff(qrs_indexes) / fs;
ibi_time = (time(qrs_indexes(1:end-1)) + ...
    time(qrs_indexes(2:end))) / 2;

% Reject intervals outside the same physiological BPM range used by the
% phase-level ECG pipeline, then remove large deviations from the median.
valid = isfinite(ibi) & ibi >= 60/160 & ibi <= 60/40;
if any(valid)
    median_ibi = median(ibi(valid));
    valid = valid & ibi >= 0.5*median_ibi & ibi <= 1.5*median_ibi;
end

ibi = ibi(valid);
ibi_time = ibi_time(valid);
if numel(ibi) < 3
    error('MultimodalPhysioKit:ECG:InsufficientGlobalIBI', ...
        'Too few valid global IBIs remain for STFT analysis.')
end

context = struct();
context.IBI = ibi;
context.IBITime = ibi_time;
context.TimeRange = [time(1) time(end)];
phase_starts = time(i_phases(1,:));
phase_ends = time(i_phases(2,:));
context.PhaseWindows = [phase_starts(:) phase_ends(:)];
end
