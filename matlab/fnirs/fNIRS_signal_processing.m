function [features_oxy,features_deoxy,intermediates] = ...
        fNIRS_signal_processing( ...
        signal_red,signal_infrared,time,fs,i_phases,age_years,gr,varargin)

current_data = fnirs_preprocess_currents( ...
    signal_red,signal_infrared,time,fs,i_phases);

fprintf('fNIRS processing in progress. Please wait..\n')
%% Initialization of flags
check_flag = gr;

% Shift time so the first sample is zero for diagnostic plots.
time_new = time-time(1);
fs_new = current_data.SamplingFrequency;
downsampling_factor = current_data.DownsamplingFactor;
downsampled_red_signal = current_data.DownsampledRed;
downsampled_infrared_signal = current_data.DownsampledInfrared;
time_downsampling = current_data.Time;
i_phases_downsampling = current_data.PhaseIndexes;
signal_filter_low_red = current_data.FilteredRed;
signal_filter_low_infrared = current_data.FilteredInfrared;

% Plot input fNIRS current signals with phase markers.
if gr
    figure
    subplot(2,1,1);
    plot(time_new,signal_red,'r','LineWidth',1.5)
    hold on
    xline(time_new(i_phases(1,:)),'--k')
    xline(time_new(i_phases(2,:)),'--k')
    xlim([0 time_new(end)])
    ylabel('fNIRS [\muA]')
    grid on
    set(gca,'FontSize',14)
    xlabel('Time (s)')
    title('Channel 1 - \lambda_1')

    subplot(2,1,2);
    plot(time_new,signal_infrared,'LineWidth',1.5)
    hold on
    xline(time_new(i_phases(1,:)),'--k')
    xline(time_new(i_phases(2,:)),'--k')
    xlim([0 time_new(end)])
    ylabel('fNIRS [\muA]')
    grid on
    set(gca,'FontSize',14)
    xlabel('Time (s)')
    title('Channel 2 - \lambda_2')
end

if check_flag
    % Plotting signals before and after downsampling with a moving average filter
    % Display a representative portion of the signal without assuming a
    % minimum recording length.
    start_sample = 1;
    end_sample = min(length(signal_infrared), max(start_sample, round(0.5*length(signal_infrared))));
    downsampled_start = max(1, round(start_sample/downsampling_factor)+1);
    downsampled_start = min(downsampled_start, length(time_downsampling));
    downsampled_end = max(downsampled_start, min(length(time_downsampling), round(end_sample/downsampling_factor)));

    figure
    plot(time_new(start_sample:end_sample),signal_infrared(start_sample:end_sample)),hold on,
    plot(time_downsampling(downsampled_start:downsampled_end),downsampled_infrared_signal(downsampled_start:downsampled_end), 'k','LineWidth',1)
    plot(time_downsampling(downsampled_start:downsampled_end),downsampled_infrared_signal(downsampled_start:downsampled_end),'o', 'MarkerFaceColor', 'black', 'MarkerSize', 4)
    %xline(time_downsampling(i_phases_downsampling(1,:)),'--k')
    %xline(time_downsampling(i_phases_downsampling(2,:)),'--k')
    legend('Original signal','Signal after the moving average filter')
    title(['Downsampling to ', num2str(fs_new), ' Hz: Infrared fNIRS signal - \lambda_2'])
    xlabel('Time [s]')
    ylabel('fNIRS [\muA]')
    axis tight
    grid on

    figure,plot(time_new(start_sample:end_sample),signal_red(start_sample:end_sample)),
    hold on, plot(time_downsampling(downsampled_start:downsampled_end),downsampled_red_signal(downsampled_start:downsampled_end), 'k','LineWidth',1)
    plot(time_downsampling(downsampled_start:downsampled_end),downsampled_red_signal(downsampled_start:downsampled_end),'o', 'MarkerFaceColor', 'black', 'MarkerSize', 4)
    %xline(time_downsampling(i_phases_downsampling(1,:)),'--k')
    %xline(time_downsampling(i_phases_downsampling(2,:)),'--k')
    legend('Original signal','Signal after the moving average filter')
    title(['Downsampling to ', num2str(fs_new), ' Hz: Red fNIRS signal - \lambda_1'])
    xlabel('Time [s]')
    ylabel('fNIRS [\muA]')
    axis tight
    grid on
end

%Plot of the downsampled signals before and after applying the low-pass filter
if gr
    figure
    plot(time_downsampling,downsampled_infrared_signal,'Color',[0 0.4470 0.7410]);
    hold on
    plot(time_downsampling,signal_filter_low_infrared,'Color',[0.6350 0.0780 0.1840]);
    xline(time_downsampling(i_phases_downsampling(1,:)),'--k')
    xline(time_downsampling(i_phases_downsampling(2,:)),'--k')
    hold off
    axis tight
    legend('Downsampled \lambda_2 signal','Signal  \lambda_2 after LPF ft=0.2Hz')
    xlabel('Time [s]')
    ylabel('fNIRS [\muA]')
    title(['Comparison before and after LPF - Infrared emitter \lambda_2 - fs = ', num2str(fs_new), ' Hz'])
    figure
    plot(time_downsampling,downsampled_red_signal,'Color',[0 0.4470 0.7410]);
    hold on
    plot(time_downsampling,signal_filter_low_red,'Color',[0.6350 0.0780 0.1840]);
    xline(time_downsampling(i_phases_downsampling(1,:)),'--k')
    xline(time_downsampling(i_phases_downsampling(2,:)),'--k')
    hold off
    axis tight
    legend('Downsampled \lambda_1 signal' ,'Signal \lambda_1 after LPF ft=0.2Hz ')
    xlabel('Time [s]')
    ylabel('fNIRS [\muA]')
    title(['Comparison before and after LPF - Red emitter \lambda_1 - fs = ', num2str(fs_new), ' Hz'])
end

% Convert filtered dual-wavelength currents to oxygenated and deoxygenated
% hemoglobin concentration changes with the modified Beer-Lambert model.
[deltaC_oxy,deltaC_deoxy,baseline] = LambertBeer( ...
    signal_filter_low_red,signal_filter_low_infrared, ...
    time_downsampling,fs_new,i_phases_downsampling,age_years,0, ...
    varargin{:});

% Plot hemoglobin concentration changes across the requested phases.
if check_flag==1 && size(i_phases,2) == 1
    for i=1:size(i_phases,2)
        figure()
        plot(time_downsampling(i_phases_downsampling(1,i):i_phases_downsampling(2,i)), deltaC_oxy(i_phases_downsampling(1,i):i_phases_downsampling(2,i)), 'Color', [0.5 0.5 0.5])
        hold on
        plot(time_downsampling(i_phases_downsampling(1,i):i_phases_downsampling(2,i)), deltaC_deoxy(i_phases_downsampling(1,i):i_phases_downsampling(2,i)),'k')
        lgd=legend('\DeltaHbO_2','\DeltaHHb');

        lgd.FontSize = 12;
        lgd.FontWeight = 'bold';
        axis tight
        hold off
        ylabel('Concentration change [mM]')
        xlabel('Time [s]')
        title('fNIRS signal Concentration change')
    end
elseif check_flag == 1 && size(i_phases,2) > 1
    for i=1:size(i_phases,2)
        figure()
        plot(time_downsampling(i_phases_downsampling(1,i):i_phases_downsampling(2,i)), deltaC_oxy(i_phases_downsampling(1,i):i_phases_downsampling(2,i)), 'Color', [0.5 0.5 0.5])
        hold on
        plot(time_downsampling(i_phases_downsampling(1,i):i_phases_downsampling(2,i)), deltaC_deoxy(i_phases_downsampling(1,i):i_phases_downsampling(2,i)),'k')
        lgd=legend('\DeltaHbO_2','\DeltaHHb');

        lgd.FontSize = 12;
        lgd.FontWeight = 'bold';
        axis tight
        hold off
        ylabel('Concentration change [mM]')
        xlabel('Time [s]')
        title('fNIRS signal Concentration change')
    end
end

% Extract features from the oxygenated hemoglobin concentration signal.
[features_oxy,oxy_diagnostics] = features_extraction( ...
    deltaC_oxy,time_downsampling,i_phases_downsampling,fs_new,0);
% Extract features from the deoxygenated hemoglobin concentration signal.
[features_deoxy,deoxy_diagnostics] = features_extraction( ...
    deltaC_deoxy,time_downsampling,i_phases_downsampling,fs_new,0);

intermediates = struct( ...
    'downsampled_red',downsampled_red_signal, ...
    'downsampled_infrared',downsampled_infrared_signal, ...
    'downsampled_time',time_downsampling, ...
    'processed_phase_indexes',i_phases_downsampling, ...
    'processed_sampling_frequency',fs_new, ...
    'filtered_red',signal_filter_low_red, ...
    'filtered_infrared',signal_filter_low_infrared, ...
    'baseline_mode',baseline.Mode, ...
    'baseline_red_value',baseline.RedValue, ...
    'baseline_infrared_value',baseline.InfraredValue, ...
    'baseline_start_index',baseline.StartIndex, ...
    'baseline_end_index',baseline.EndIndex, ...
    'oxy_concentration',deltaC_oxy, ...
    'deoxy_concentration',deltaC_deoxy, ...
    'filtered_oxy_concentration',oxy_diagnostics.FilteredSignal, ...
    'filtered_deoxy_concentration',deoxy_diagnostics.FilteredSignal);
end
