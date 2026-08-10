function temperature_features = temperature_signal_processing(temp,time,i_phases,fs,gr)
%TEMPERATURE_SIGNAL_PROCESSING Extract phase-level skin-temperature features.
%
% Inputs
% ------
% TEMP     : Skin-temperature signal in degrees Celsius.
% TIME     : Time vector in seconds, same length as TEMP.
% I_PHASES : 2-by-N matrix of one-based inclusive phase limits.
% FS       : Sampling frequency in Hz.
% GR       : Logical flag. If true, diagnostic figures are generated.
%
% Output
% ------
% TEMPERATURE_FEATURES : Cell array containing feature names and phase-level
% values.

fprintf('Temperature processing in progress. Please wait..\n')

%% Temperature filtering process
temp_f = filtering(temp,fs,0,1,0,5,0);

%% Signal processing
temperature_features = {
    'temperature_delta';
    'temperature_mean';
    'temperature_std';
    'temperature_delta_over_time';
    'temperature_linear_slope';
    'temperature_derivative_mean';
    'temperature_derivative_std';
    'temperature_derivative_linear_slope'
    };
% Start signal processing
all_indexes_phases = split_signal_phases(i_phases);

for phase = 1 : length(all_indexes_phases)
    max_rel_tot = [];
    min_rel_tot = [];
    indexes_phase = all_indexes_phases{phase};
    filtered_temperature_phase = temp_f(indexes_phase);
    diff_filtered_temperature_phase = diff(filtered_temperature_phase);
    time_phase = time(indexes_phase);
    mean_filtered_temperature_phase = mean(filtered_temperature_phase);
    mean_diff_filtered_temperature_phase = mean(diff_filtered_temperature_phase);
    % Screen short windows for abrupt artifacts before feature extraction.
    window_length = 2*fs;                                                   % in samples (-)
    threshold_percentage_max_relative = 1.06;                               % percentage (%)
    threshold_diff = 100;                                                   % percentage (%)
    length_filtered_temperature_phase = length(filtered_temperature_phase); % in samples
    if window_length < length_filtered_temperature_phase
        indexes_windows_phase = 1:window_length:length_filtered_temperature_phase;
        if length_filtered_temperature_phase - indexes_windows_phase(end) < window_length && ...
                length_filtered_temperature_phase - indexes_windows_phase(end) ~= 0
            indexes_windows_phase(end+1) = length_filtered_temperature_phase;
        end
    else
        indexes_windows_phase = [1 length_filtered_temperature_phase];
    end
    for j = 1:length(indexes_windows_phase)-1
        interval_samples = indexes_windows_phase(j):indexes_windows_phase(j+1);
        filtered_temperature_phase_window = filtered_temperature_phase(interval_samples);
        diff_filtered_temperature_phase_window = diff(filtered_temperature_phase_window);
        mean_diff_filtered_temperature_phase_window = mean(abs(diff_filtered_temperature_phase_window));
        % Analysis with relative peaks detection
        [indexes_max_rel,indexes_min_rel] = find_max_min(filtered_temperature_phase_window,0,0);
        interval_indexes_max_rel = interval_samples(indexes_max_rel);
        interval_indexes_min_rel = interval_samples(indexes_min_rel);
        max_rel_tot = [max_rel_tot interval_indexes_max_rel];
        min_rel_tot = [min_rel_tot interval_indexes_min_rel];
        max_min_indexes = [indexes_max_rel indexes_min_rel];
        if isempty(filtered_temperature_phase_window(max_min_indexes))
            max_min_indexes = [1 round(length(interval_samples)/2) length(interval_samples)];
        end
        if (any(filtered_temperature_phase_window(max_min_indexes) > threshold_percentage_max_relative * mean_filtered_temperature_phase) || ...
                any(filtered_temperature_phase_window(max_min_indexes) < (2-threshold_percentage_max_relative) * mean_filtered_temperature_phase)) && ...
                mean_diff_filtered_temperature_phase_window > threshold_diff * mean_diff_filtered_temperature_phase
            filtered_temperature_phase(interval_samples,2) = zeros(length(interval_samples),1);
        else
            filtered_temperature_phase(interval_samples,2) = ones(length(interval_samples),1);
        end
    end

    % Variables update
    time_phase_final = time_phase(filtered_temperature_phase(:,2) == 1);
    filtered_temperature_phase_final = filtered_temperature_phase(filtered_temperature_phase(:,2) == 1,1);
    diff_filtered_temperature_phase_final = diff(filtered_temperature_phase_final);
    diff_time_phase_final = diff(time_phase_final);
    derivative_temperature_phase_final = diff_filtered_temperature_phase_final./diff_time_phase_final;

    %% Features extraction section
    % Delta temperature
    temperature_features{1,phase+1} = filtered_temperature_phase_final(end) - filtered_temperature_phase_final(1);

    % Mean temperature
    temperature_features{2,phase+1} = mean(filtered_temperature_phase_final);

    % Standard deviation of temperature
    temperature_features{3,phase+1} = std(filtered_temperature_phase_final);

    % Variation of the temperature over time
    temperature_features{4,phase+1} = temperature_features{1,phase+1}/(time_phase_final(end)-time_phase_final(1));

    % Linear slope of the temperature in the phase.
    c1 = polyfit(time_phase_final,filtered_temperature_phase_final,1);
    temperature_features{5,phase+1} = c1(1);

    % Mean temperature first derivative
    temperature_features{6,phase+1} = mean(derivative_temperature_phase_final);

    % Standard deviation of temperature first derivative
    temperature_features{7,phase+1} = std(derivative_temperature_phase_final);

    % Interpolated slope of first derivative
    c2 = polyfit(time_phase_final(2:end),derivative_temperature_phase_final,1);
    temperature_features{8,phase+1} = c2(1);

end

if gr
    % Plot Temperature vs Time
    figure('Name','Temperature vs time')
    hold on
    plot(time,temp,'LineWidth',1.5)
    plot(time,temp_f,'LineWidth',1.5)
    xline(time([i_phases(1,:),i_phases(2,:)]))
    grid on
    xlabel('Time (s)')
    ylabel('Raw temperature (°C)')
    set(gca,'FontSize',15)
    % Plot First Derivative vs Time
    sign_diff_filtered_temperature_phase = sign(diff(temp)) >= 0;
    figure('Name','Temperature First Derivative vs Time')
    hold on
    plot(time,temp,'LineWidth',1.5)
    derivative_time = time(2:end);
    derivative_temperature = temp(2:end);
    plot(derivative_time(sign_diff_filtered_temperature_phase == 0),derivative_temperature(sign_diff_filtered_temperature_phase == 0),'*r','LineWidth',0.5)
    plot(derivative_time(sign_diff_filtered_temperature_phase == 1),derivative_temperature(sign_diff_filtered_temperature_phase == 1),'*b','LineWidth',0.5)
    xline(time([i_phases(1,:),i_phases(2,:)]))
    grid on
    xlabel('Time (s)')
    ylabel('Raw temperature (°C)')
    set(gca,'FontSize',15)
end




end
