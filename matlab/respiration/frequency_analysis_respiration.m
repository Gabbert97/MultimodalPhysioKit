function [mean_PLF_final,std_PLF_final,mean_PHF_final,std_PHF_final,mean_LF_HF_ratio_final,std_LF_HF_ratio_final] = frequency_analysis_respiration(respiration_rate,length_rr_phase,gr)
%FREQUENCY_ANALYSIS_RESPIRATION Estimate phase-level respiration-rate spectra.
%
% Inputs
% ------
% RESPIRATION_RATE : Breath-rate sequence in breaths per minute.
% LENGTH_RR_PHASE  : Number of respiration-rate samples assigned to each phase.
% GR               : Logical flag. If true, diagnostic figures are generated.
%
% Outputs
% -------
% Phase-level low-frequency, high-frequency, and LF/HF summary values.

%% Frequency-analysis configuration
% 1 = parametric analysis; 2 = STFT.
frequency_analysis_label = 1;
% Parametric PSD method: 1 = Burg, 2 = modified covariance,
% 3 = Yule-Walker, 4 = Welch.
label_PSD = 1;
% Moving-window length in respiration-rate samples.
winlen_interval = 200;
final_window = winlen_interval(end);

%% Parametric frequency-analysis parameters
N = 50;                                                                     % Max order for parametric evaluations
perc=5;
NFFT = 512;
% Enable internal order-selection diagnostics during development.
check_flag = 0;

%% STFT frequency-analysis parameters
nfft = 256;

%% Index initialization process
indexes_end_phase = cumsum(length_rr_phase);
indexes_end_phase = [0;indexes_end_phase];
% Convert breath rate to breath interval in seconds.
signal = (60./respiration_rate)';      % [s]
f_mean = 1/mean(signal);            % [1/s]

%% Time-Frequency Analysis
for j = 1:length(winlen_interval)
    winlen = winlen_interval(j);
    %% Parametric analysis
    if frequency_analysis_label == 1
        % Zero-pad the signal so each output sample has a centered window.
        if size(respiration_rate,2)==1
            respiration_rate = respiration_rate.';
        end
        if rem(winlen,2) == 0
            zero_padding_signal=[zeros(1,winlen/2-1) signal zeros(1,winlen/2)];
        else
            zero_padding_signal=[zeros(1,(winlen-1)/2) signal zeros(1,(winlen-1)/2)];
        end
        for i = 1 : length(signal)
            % Windowed breath-interval segment.
            windowed_signal = zero_padding_signal(i:winlen-1+i);
            % Estimate the nonzero window mean introduced by edge padding.
            zeros_number = length(windowed_signal(windowed_signal == 0));
            mean_windowed_signal = sum(windowed_signal)/(length(windowed_signal)-zeros_number);
            mean_frequency_windowed_signal = 1/mean_windowed_signal;
            if length(windowed_signal) <= N
                N = length(windowed_signal) - 1;
            end
            % Burg PSD method implementation
            if label_PSD == 1
                % Select pburg model order.
                for NN = 2:N
                    [~,e1(NN)] = arburg(windowed_signal,NN);
                end
                NN = 2:N;
                if check_flag == 1
                    figure
                    plot(NN,e1(2:N),'r')
                    xlabel ('order'), ylabel('variance'), axis tight
                end
                % Select the first order near the asymptotic error variance.
                asymptotic_threshold = perc*e1(end)/100 + e1(end);
                % Exclude the first variance value, which is unstable for
                % very low model orders.
                candidate_orders = find(e1(2:N)<asymptotic_threshold);
                order1 = candidate_orders(1);
                % Parametric spectral density analysis
                windowed_signal(windowed_signal ~= 0) = windowed_signal(windowed_signal ~= 0) - mean_windowed_signal;
                [Pb(:,i), w(:,i)] = pburg(windowed_signal,order1,NFFT);
                % Frequency normalization process
                f(:,i) = w(:,i)/(2*pi)*mean_frequency_windowed_signal;

                % Modified covariance PSD method.
            elseif label_PSD == 2
                for NN = 2:N
                    [~,e(NN)] = armcov(windowed_signal,NN);
                end
                if check_flag == 1
                    figure(4)
                    NN = 2:N;
                    plot(NN,e(2:N),'r')
                    xlabel ('order'), ylabel('variance'), axis tight
                end

                asymptotic_threshold = perc*e(end)/100+e(end);
                candidate_orders = find(e(2:N)<asymptotic_threshold);
                order = candidate_orders(1);

                [Pb(:,i),w(:,i)] = pmcov(windowed_signal-mean_windowed_signal,order,NFFT);
                f(:,i) = w(:,i)/(2*pi)*mean_frequency_windowed_signal;

                % Yule-Walker
            elseif label_PSD == 3
                for NN = 2:N
                    [~,e(NN)] = aryule(windowed_signal,NN);
                end
                if check_flag == 1
                    figure
                    NN = 2:N;
                    plot(NN,e(2:N),'r')
                    xlabel ('order'), ylabel('variance'), axis tight
                end
                asymptotic_threshold = perc*e(end)/100+e(end);
                candidate_orders = find(e(2:N)<asymptotic_threshold);
                order = candidate_orders(1);
                [Pb(:,i),w(:,i)] = pyulear(windowed_signal-mean_windowed_signal,order,NFFT);
                f(:,i) = w(:,i)/(2*pi)*mean_frequency_windowed_signal;
                % Non-parametric Welch PSD evaluation method
            elseif label_PSD == 4
                [Pb(:,i),w(:,i)] = pwelch(windowed_signal-mean_windowed_signal);
                f(:,i) = w(:,i)/(2*pi)*mean_frequency_windowed_signal;

            else
                error('MultimodalPhysioKit:InvalidSpectralMethod', ...
                    'Select Burg, modified covariance, Yule-Walker, or Welch spectral analysis.')
            end
        end

        %% STFT Frequency Analysis
    elseif frequency_analysis_label == 2
        %% Short-time Fourier transform analysis
        [S1,~,f] = STFT(signal,nfft,f_mean,winlen(j),hamming(winlen(j)));
        Pb = abs(S1);
    end
end
length_time = length(signal);
t_plot = 1:length_time;
% Power spectral density representation.
if gr
    figure
    hold on
    mesh(t_plot,f,abs(Pb),'FaceColor','interp')
    colormap jet
    xline(indexes_end_phase,'w--','LineWidth',0.8)
    yl1 = yline(0.15,'w','LF');
    yl1.LabelVerticalAlignment = 'bottom';
    yl1.LabelHorizontalAlignment = "left";
    yl2 = yline(0.15,'w','HF');
    yl2.LabelVerticalAlignment = 'top';
    yl2.LabelHorizontalAlignment = "left";
    xlabel('Time (samples)')
    ylabel('Frequency (Hz)')
    zlabel('PSD')
    xlim([t_plot(1) t_plot(end)])
    ylim([min(f(:)) max(f(:))])
    set(gca,'FontSize',15)
    grid on
    view(2)
    title(sprintf('Window length: %d', winlen(j)))
end
%% LF over HF evaluation
for i = 1:length_time
    index = f>=0.04 & f<=0.18;
    Ptot = sum(Pb(:,i));
    P = sum(Pb(index,i));
    PLF(j,i) = P/Ptot;
    index = f>0.18 & f<=0.4;
    P = sum(Pb(index,i));
    PHF(j,i) = P/Ptot;
    LF_HF_ratio(j,i) = PLF(j,i)/PHF(j,i);
end
% Evaluate phase-level LF, HF, and LF/HF summary statistics.
for i = 1:length(indexes_end_phase)-1
    indexes_phases = indexes_end_phase(i)+1:indexes_end_phase(i+1);
    mean_PLF(j,i) = mean(PLF(j,indexes_phases));
    std_PLF(j,i) = std(PLF(j,indexes_phases));
    mean_PHF(j,i) = mean(PHF(j,indexes_phases));
    std_PHF(j,i) = std(PHF(j,indexes_phases));
    mean_LF_HF_ratio(j,i) = mean(LF_HF_ratio(j,indexes_phases));
    std_LF_HF_ratio(j,i) = std(LF_HF_ratio(j,indexes_phases));
end
% Plot the LF/HF trend.
if gr
    figure
    hold on
    plot(1:length_time,LF_HF_ratio(j,:),'LineWidth',1.5)
    xline(indexes_end_phase,'k--','LineWidth',0.8)
    xlabel('samples (-)')
    ylabel('LF/HF (-)')
    set(gca,'FontSize',15)
end

%% Final window decision and output creation
final_window_index = winlen == final_window;
mean_PLF_final = mean_PLF(final_window_index,:);
std_PLF_final = std_PLF(final_window_index,:);
mean_PHF_final = mean_PHF(final_window_index,:);
std_PHF_final = std_PHF(final_window_index,:);
mean_LF_HF_ratio_final = mean_LF_HF_ratio(final_window_index,:);
std_LF_HF_ratio_final = std_LF_HF_ratio(final_window_index,:);

end
