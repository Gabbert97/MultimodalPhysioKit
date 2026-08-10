function [signal_downsampling,time_downsampling, i_phases_downsampling] = downsampling(signal,time,fs,fs_new,i_phases, gr)
%DOWNSAMPLING Downsample a signal with non-overlapping moving averages.
%
% Inputs
% ------
% SIGNAL   : Input signal vector.
% TIME     : Time vector in seconds, same length as SIGNAL.
% FS       : Original sampling frequency in Hz.
% FS_NEW   : Target sampling frequency in Hz.
% I_PHASES : 2-by-N matrix of phase limits in samples.
% GR       : Logical flag. If true, diagnostic figures are generated.

% Calculate the downsampling factor.
downsampling_factor = fs/ fs_new;
N = downsampling_factor;

%% Downsampling phases
% Convert phase indexes to the downsampled sample grid.
i_phases_downsampling = round(i_phases./downsampling_factor);
%% Downsampling signal
% Initialize downsampled signal and time vectors.
signal_downsampling = zeros(round(length(signal)/downsampling_factor),1);
time_downsampling = zeros(length(signal_downsampling),1);

% Average non-overlapping groups of samples and assign each value to the
% central time point of the group.
for i = 1:length(signal_downsampling)
    start_index = (i-1)*downsampling_factor + 1;
    end_index = start_index + downsampling_factor-1;

    % Check if the end index exceeds the length of the original signal
    if end_index>length(signal)    
        % Adjust the end index for the last group
        end_index = length(signal);
    end
    x_grouped = signal(start_index:end_index);
    signal_downsampling(i,1) = mean(x_grouped);
    time_index = start_index + round(length(x_grouped)/2) - 1;
    time_downsampling(i,1) = time(time_index);
end

% Plot the original and downsampled signals.
if gr==1
    figure;
    subplot(2,1,1);
    plot(time,signal);
    xlabel('Time (s)');
    ylabel('fNIRS [\muA]')
    title(sprintf('Original signal sampled at %.3g Hz', fs))
    xline(time(i_phases(1,:)),'--')
    xline(time(i_phases(2,:)),'--')
    ylim([min(signal) max(signal)])
    xlim([0 time(end)])
    
    subplot(2,1,2);
    plot(time_downsampling, signal_downsampling);
    xline(time_downsampling(i_phases_downsampling(1,:)),'--')
    xline(time_downsampling(i_phases_downsampling(2,:)),'--')
    xlabel('Time (s)');
    ylabel('fNIRS [\muA]')
    title(sprintf('Downsampled signal at %.3g Hz', fs_new))
    ylim([min(signal) max(signal)])
    xlim([0 time_downsampling(end)])

    figure;
    plot(time,signal);
    hold on
    plot(time_downsampling, signal_downsampling, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    xline(time_downsampling(i_phases_downsampling(1,:)),'--')
    xline(time_downsampling(i_phases_downsampling(2,:)),'--')
    hold off
    xlabel('Time (s)');
    ylabel('fNIRS [\muA]')
    ylim([min(signal) max(signal)])
    xlim([0 time(end)])
    legend('Original signal','Downsampled signal')
    title('Comparison between the signal before and after downsampling with a moving average filter')

end
if gr==1 && length(signal_downsampling) > 1
    start_sample = 1;
    end_sample = min(length(signal_downsampling), ...
        max(start_sample + 1, round(0.5*length(signal_downsampling))));
    original_start = max(1, start_sample*N);
    original_end = min(length(signal), end_sample*N);

    figure
    plot(time(original_start:original_end),signal(original_start:original_end),'b');
    hold on
    plot(time_downsampling(start_sample:end_sample),signal_downsampling(start_sample:end_sample),'r')
    title(sprintf('Downsampling at %.3g Hz using mean values', fs_new))
end
end
