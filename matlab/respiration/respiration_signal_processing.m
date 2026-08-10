function respiration_features = respiration_signal_processing(respiration,time,i_phases,fs,gr)
if mod(fs,100) ~= 0
    error('MultimodalPhysioKit:InvalidSamplingFrequency', ...
        'Respiration sampling frequency must be an integer multiple of 100 Hz.')
end

fprintf('Respiration processing in progress. Please wait..\n')
%% Variables initialization
fs_new = 100;                                                               % Downsample new frequency [Hz]
i_min_def = [];
i_max_def = [];                                                             % Maximum points array for wave detection
respiration_rate_tot = [];                                                  % Respiration rate tot
length_respiration_rate_tot = [];                                           % Length respiration rate
respiration_features = {'respiration_breath_rate_mean';'respiration_breath_rate_std';'respiration_inspiratory_time_mean';'respiration_inspiratory_time_std';'respiration_expiratory_time_mean';...
    'respiration_expiratory_time_std';'respiration_timing_ratio_mean';'respiration_timing_ratio_std';'respiration_amplitude_mean';'respiration_amplitude_std';...
    'respiration_minute_ventilation_mean';'respiration_minute_ventilation_std';'respiration_coefficient_of_variation';'respiration_autocorrelation'};

% Check flag
check_flag = 0;                                                             % [1] - plot on ; [0] - plot off

%% Downsample parameters
ts_new = 1/fs_new;                                                                  % Sampling time [s]
down = fs/fs_new;
respiration_downsampled = downsample(respiration,down);
time_downsampled = downsample(time,down);
i_phases_downsampled = ceil(i_phases./down);

%% Indexes phase evaluation
all_indexes_phases = split_signal_phases(i_phases_downsampled);

%% Diagnostic plot
if gr
    figure
    plot(time_downsampled,respiration_downsampled,'LineWidth',1.5)
    hold on
    xline(time_downsampled(i_phases_downsampled(1,:)),'--')
    xline(time_downsampled(i_phases_downsampled(2,:)),'--')
    xlabel('Time [s]')
    ylabel('Respiration [s]')
    set(gca,'FontSize',14)
end
%% Analysis part
% respiration_downsampled = respiration_downsampled - mean(respiration_downsampled);
for phase = 1 : length(all_indexes_phases)
    % Phase indexes definition
    indexes_phase = all_indexes_phases{phase};
    % The respiration signal is now taken into account and analysed in the
    % considered i-th phase
    respiration_downsampled_phase = respiration_downsampled(indexes_phase);
    time_downsampled_phase = time_downsampled(indexes_phase);
    filtered_phase_respiration = filtering(respiration_downsampled_phase,fs_new,0,1.1,0,4,0);
    %% Window loop to extract all the characteristics from the respiration
    % It was decided to work on the filtered signal in the phase because it
    % retains the respiration-band frequency content.
    [i_filt_peaks_min,i_filt_peaks_max] = find_peaks_respiration(filtered_phase_respiration,0,1);
    % Control start-end minimum and maximum points
    if i_filt_peaks_min(1) >= i_filt_peaks_max(1)
        i_filt_peaks_max = i_filt_peaks_max(2:end);
    end
    if i_filt_peaks_min(end) <= i_filt_peaks_max(end)
        i_filt_peaks_max = i_filt_peaks_max(1:end-1);
    end
    % Check plot
    if gr
        hold on
        plot(time_downsampled_phase,filtered_phase_respiration,'LineWidth',1.5)
        plot(time_downsampled_phase(i_filt_peaks_min),filtered_phase_respiration(i_filt_peaks_min),'o','LineWidth',1.5)
        plot(time_downsampled_phase(i_filt_peaks_max),filtered_phase_respiration(i_filt_peaks_max),'*','LineWidth',1.5)
    end
    %% Control section
    count_distance_def = 1;                                                 % Counter distances definitive [-]
    minimum_time_distance_def = [];
    max_min_distance_down_def = [];
    max_min_distance_up_def = [];
    i_min_def_phase = [];
    i_max_def_phase = [];                                                             % Maximum points array for wave detection

    % Minimum points time distance
    minimum_time_distance = diff(time_downsampled_phase(i_filt_peaks_min));
    mean_time_distance = mean(minimum_time_distance);
    % Minimum points amplitude position
    minimum_amplitude_distance = abs(diff(filtered_phase_respiration(i_filt_peaks_min)));
    mean_amplitude_distance = mean(minimum_amplitude_distance);
    % Max - min distances up & down
    max_min_distance_down = abs(filtered_phase_respiration(i_filt_peaks_min(1:end-1))-filtered_phase_respiration(i_filt_peaks_max));
    max_min_distance_up = abs(filtered_phase_respiration(i_filt_peaks_min(2:end))-filtered_phase_respiration(i_filt_peaks_max));
    mean_max_min_distance = mean([max_min_distance_down;max_min_distance_up]);
    mean_max_min_distance_diff = mean(abs(max_min_distance_up-max_min_distance_down));

    % Starting point definition
    i_start_minimum_point = find(max_min_distance_down >= 0.6 * mean_max_min_distance);
    i_min_def_phase(1) = i_filt_peaks_min(i_start_minimum_point(1));

    % Check of the intervals
    for i = i_start_minimum_point(1) : length(minimum_time_distance)
        % Definiton of the wave parameters -> 2 minimum point and 1
        % relative max
        i_min_point_down = i_min_def_phase(count_distance_def);
        i_min_point_up = i_filt_peaks_min(i+1);
        % max point evaluation
        [~,i_max_relative_point] = max(filtered_phase_respiration(i_min_point_down:i_min_point_up));
        i_max_point = i_min_point_down + i_max_relative_point - 1;
        if check_flag
            figure
            hold on
            plot(1:length(filtered_phase_respiration),filtered_phase_respiration,'b','Linewidth',1.5)
            plot(i_min_point_down:i_min_point_up,filtered_phase_respiration(i_min_point_down:i_min_point_up),'r','LineWidth',1.5)
            plot([i_min_point_down,i_min_point_up,],filtered_phase_respiration([i_min_point_down,i_min_point_up,]),'o')
            plot(i_max_point,filtered_phase_respiration(i_max_point),'*')
        end

        if i < length(minimum_time_distance)
            % Wave parameters evaluation
            minimum_time_distance_wave = time_downsampled_phase(i_min_point_up) - time_downsampled_phase(i_min_point_down);
            %             minimum_time_distance_wave_wave = time_downsampled_phase(i_filt_peaks_min(i+2)) - time_downsampled_phase(i_min_point_up);
            max_min_distance_down_wave = filtered_phase_respiration(i_max_point) - filtered_phase_respiration(i_min_point_down);
            max_min_distance_up_wave = filtered_phase_respiration(i_max_point) - filtered_phase_respiration(i_min_point_up);
            max_min_distance_up_wave_wave = abs(filtered_phase_respiration(i_filt_peaks_max(i+1))-filtered_phase_respiration(i_min_point_up));

            % Wave control and definitive minimmum points evaluation
            if minimum_time_distance_wave >= 0.1 * mean_time_distance && max_min_distance_down_wave > 0.2 * mean_max_min_distance && ...
                    max_min_distance_up_wave > 0.2 * mean_max_min_distance && max_min_distance_up_wave_wave > 0.2 * mean_max_min_distance

                i_min_def_phase = [i_min_def_phase;i_min_point_up];
                i_max_def_phase = [i_max_def_phase;i_max_point];
                minimum_time_distance_def = [minimum_time_distance_def;minimum_time_distance_wave];
                max_min_distance_down_def = [max_min_distance_down_def;max_min_distance_down_wave];
                max_min_distance_up_def = [max_min_distance_up_def;max_min_distance_up_wave];

                if length(i_min_def_phase) > 7
                    mean_interval = (count_distance_def-6):count_distance_def;
                    mean_time_distance = mean(minimum_time_distance_def(mean_interval));
                    mean_max_min_distance = mean([max_min_distance_down_def(mean_interval);max_min_distance_up_def(mean_interval)]);
                end

                count_distance_def = count_distance_def + 1;
            end

        else
            % Wave parameters evaluation
            minimum_time_distance_wave = time_downsampled_phase(i_min_point_up) - time_downsampled_phase(i_min_point_down);
            max_min_distance_down_wave = filtered_phase_respiration(i_max_point) - filtered_phase_respiration(i_min_point_down);
            max_min_distance_up_wave = filtered_phase_respiration(i_max_point) - filtered_phase_respiration(i_min_point_up);


            if minimum_time_distance_wave >= 0.5 * mean_time_distance && max_min_distance_down_wave > 0.4 * mean_max_min_distance && ...
                    max_min_distance_up_wave > 0.4 * mean_max_min_distance

                i_min_def_phase = [i_min_def_phase;i_min_point_up];
                i_max_def_phase = [i_max_def_phase;i_max_point];
                minimum_time_distance_def = [minimum_time_distance_def;minimum_time_distance_wave];
                max_min_distance_down_def = [max_min_distance_down_def;max_min_distance_down_wave];
                max_min_distance_up_def = [max_min_distance_up_def;max_min_distance_up_wave];

                if length(i_min_def_phase) > 2
                    mean_interval = (count_distance_def-1):count_distance_def;
                    mean_time_distance = mean(minimum_time_distance_def(mean_interval));
                    mean_max_min_distance = mean([max_min_distance_down_def(mean_interval);max_min_distance_up_def(mean_interval)]);
                end

                count_distance_def = count_distance_def + 1;
            end
        end
    end
    if check_flag
        figure
        hold on
        plot(time_downsampled_phase,filtered_phase_respiration,'LineWidth',1.5)
        plot(time_downsampled_phase(i_min_def_phase),filtered_phase_respiration(i_min_def_phase),'o')
        plot(time_downsampled_phase(i_max_def_phase),filtered_phase_respiration(i_max_def_phase),'*')
        grid on
        xlabel('Time [s]')
        ylabel('Chest Dispacement [%]')
    end

    i_min_def = [i_min_def;i_min_def_phase+indexes_phase(1)-1];
    i_max_def = [i_max_def;i_max_def_phase+indexes_phase(1)-1];             % Maximum points array for wave detection

    %% Features extraction and evaluation
    % Features extraction
    diff_t_max = diff(time_downsampled_phase(i_max_def_phase));                % Differences between maximum points in time intervals [s]
    respiration_rate = 60./diff_t_max;
    % Control
    indexes_rr_ok = find(respiration_rate >= 5 & respiration_rate <= 30);
    i_max_rr_ok = [i_max_def_phase(1);i_max_def_phase(indexes_rr_ok+1)];
    i_min_rr_down_ok = [i_min_def_phase(1);i_min_def_phase(indexes_rr_ok+1)];
    i_min_rr_up_ok = [i_min_def_phase(2);i_min_def_phase(indexes_rr_ok+2)];
    length_wave_ok = length(i_max_rr_ok);
    if check_flag
        figure
        hold on
        plot(time_downsampled_phase,filtered_phase_respiration,'LineWidth',1.5)
        plot(time_downsampled_phase(i_min_rr_down_ok),filtered_phase_respiration(i_min_rr_down_ok),'b*','LineWidth',2)
        plot(time_downsampled_phase(i_min_rr_up_ok),filtered_phase_respiration(i_min_rr_up_ok),'ro')
        plot(time_downsampled_phase(i_max_rr_ok),filtered_phase_respiration(i_max_rr_ok),'y*')
        grid on
        xlabel('Time [s]')
        ylabel('Chest Dispacement [%]')
    end
    % Retain rates within the configured plausible screening range. This
    % numerical range is a project processing choice and should be reviewed
    % for the population and protocol under study.
    %% Respiration rate
    respiration_rate = respiration_rate(indexes_rr_ok);
    if check_flag
        figure
        plot(1:length(respiration_rate),respiration_rate,'LineWidth',1.5)
        grid on
        xlabel('Samples [-]')
        ylabel('Breath Rate [-]')
    end
    % Respiration rate mean
    respiration_features{1,phase+1} = mean(respiration_rate);
    % Respiration rate standard deviation
    respiration_features{2,phase+1} = std(respiration_rate);
    %% Inspiratory & expiratory time & timing ratio
    inspiratory_time = time_downsampled_phase(i_max_rr_ok) - time_downsampled_phase(i_min_rr_down_ok);
    expiratory_time = time_downsampled_phase(i_min_rr_up_ok) - time_downsampled_phase(i_max_rr_ok);
    respiration_features{3,phase+1} = mean(inspiratory_time);
    respiration_features{4,phase+1} = std(inspiratory_time);
    respiration_features{5,phase+1} = mean(expiratory_time);
    respiration_features{6,phase+1} = std(expiratory_time);
    respiration_features{7,phase+1} = mean(inspiratory_time./expiratory_time);
    respiration_features{8,phase+1} = std(inspiratory_time./expiratory_time);
    %% Tidal Volume
    lower_point_amplitude = [];
    for j = 1:length_wave_ok
        minimum_time_point_wave_ok = [time_downsampled_phase(i_min_rr_down_ok(j));time_downsampled_phase(i_min_rr_up_ok(j))];
        minimum_filtered_respiration_wave_ok = [filtered_phase_respiration(i_min_rr_down_ok(j));filtered_phase_respiration(i_min_rr_up_ok(j))];
        polinomial_coefficient = polyfit(minimum_time_point_wave_ok,minimum_filtered_respiration_wave_ok,1);
        lower_point_amplitude(j) = polyval(polinomial_coefficient,time_downsampled_phase(i_max_rr_ok(j)));
    end
    tidal_volume = filtered_phase_respiration(i_max_rr_ok)-lower_point_amplitude';
    respiration_features{9,phase+1} = mean(tidal_volume);
    respiration_features{10,phase+1} = std(tidal_volume);
    %% Minute Ventilation
    minute_ventilation = respiration_rate.*tidal_volume(2:end);
    respiration_features{11,phase+1} = mean(minute_ventilation);
    respiration_features{12,phase+1} = std(minute_ventilation);
    %% Respiration Rate Variability
    % Respiration coefficient of variation
    respiration_features{13,phase+1} = (std(respiration_rate)/mean(respiration_rate))*100;
    % Autocorrelation of the respiration rate
    respiration_rate_1 = respiration_rate(1:end-1);
    respiration_rate_2 = respiration_rate(2:end);
    respiration_features{14,phase+1} = abs(corr(respiration_rate_1, respiration_rate_2));
    %% Upload of the overall respiration rate array
    % It is necessary for the following frequency analysis
    respiration_rate_tot = [respiration_rate_tot;respiration_rate];
    length_respiration_rate_tot = [length_respiration_rate_tot;length(respiration_rate)];
end
%% Frequency Analysis
%[mean_PLF_final,std_PLF_final,mean_PHF_final,std_PHF_final,mean_LF_HF_ratio_final,std_LF_HF_ratio_final] = frequency_analysis_respiration(respiration_rate_tot,length_respiration_rate_tot,gr);
%for phase = 1 : length(all_indexes_phases)
%    respiration_features{13,phase+1} = mean_PLF_final(phase);
%     respiration_features{14,phase+1} = std_PLF_final(phase);
%     respiration_features{15,phase+1} = mean_PHF_final(phase);
%     respiration_features{16,phase+1} = std_PHF_final(phase);
%     respiration_features{17,phase+1} = mean_LF_HF_ratio_final(phase);
%     respiration_features{18,phase+1} = std_LF_HF_ratio_final(phase);
% end
end
