function [PLF,PHF,LF_HF_ratio] = parametric_frequency_analysis(signal,winlen,indexes_end_phase,label_PSD,output_label,gr)
%PARAMETRIC_FREQUENCY_ANALYSIS Estimate moving LF and HF spectral power.

check_flag = 0;
saving_flag = 0;
figure_output_folder = 'figures';
if saving_flag && ~isfolder(figure_output_folder)
    mkdir(figure_output_folder)
end

if size(signal,2) == 1
    signal = signal.';
end
if rem(winlen,2) == 0
    zero_padding_signal = [zeros(1,winlen/2-1) signal zeros(1,winlen/2)];
else
    zero_padding_signal = [zeros(1,(winlen-1)/2) signal zeros(1,(winlen-1)/2)];
end

N = 50;
perc = 5;
NFFT = 512;
for i = 1:length(signal)
    windowed_signal = zero_padding_signal(i:winlen-1+i);
    zeros_number = length(windowed_signal(windowed_signal == 0));
    mean_windowed_signal = sum(windowed_signal)/(length(windowed_signal)-zeros_number);
    mean_frequency_windowed_signal = 1/mean_windowed_signal;
    if length(windowed_signal) <= N
        N = length(windowed_signal) - 1;
    end

    if label_PSD == 1
        for NN = 2:N
            [~,e1(NN)] = arburg(windowed_signal,NN);
        end
        NN = 2:N;
        if check_flag == 1
            figure
            plot(NN,e1(2:N),'r')
            xlabel('Order')
            ylabel('Variance')
            axis tight
        end
        asymptotic_threshold = perc*e1(end)/100 + e1(end);
        candidate_orders = find(e1(2:N) < asymptotic_threshold);
        order = candidate_orders(1);
        windowed_signal(windowed_signal ~= 0) = ...
            windowed_signal(windowed_signal ~= 0) - mean_windowed_signal;
        [Pb(:,i),w(:,i)] = pburg(windowed_signal,order,NFFT);
    elseif label_PSD == 2
        for NN = 2:N
            [~,e(NN)] = armcov(windowed_signal,NN);
        end
        asymptotic_threshold = perc*e(end)/100 + e(end);
        candidate_orders = find(e(2:N) < asymptotic_threshold);
        order = candidate_orders(1);
        [Pb(:,i),w(:,i)] = pmcov(windowed_signal-mean_windowed_signal,order,NFFT);
    elseif label_PSD == 3
        for NN = 2:N
            [~,e(NN)] = aryule(windowed_signal,NN);
        end
        asymptotic_threshold = perc*e(end)/100 + e(end);
        candidate_orders = find(e(2:N) < asymptotic_threshold);
        order = candidate_orders(1);
        [Pb(:,i),w(:,i)] = pyulear(windowed_signal-mean_windowed_signal,order,NFFT);
    elseif label_PSD == 4
        [Pb(:,i),w(:,i)] = pwelch(windowed_signal-mean_windowed_signal);
    else
        error('MultimodalPhysioKit:InvalidSpectralMethod', ...
            'Select Burg, modified covariance, Yule-Walker, or Welch spectral analysis.')
    end
    f(:,i) = w(:,i)/(2*pi)*mean_frequency_windowed_signal;
end

index_LF = f >= 0.04 & f <= 0.15;
index_HF = f > 0.15 & f <= 0.4;
for i = 1:length(signal)
    total_power = sum(Pb(:,i));
    PLF(i) = sum(Pb(index_LF(:,i),i))/total_power;
    PHF(i) = sum(Pb(index_HF(:,i),i))/total_power;
    LF_HF_ratio(i) = PLF(i)/PHF(i);
end

if gr
    t_plot = 1:length(signal);
    low_indexes = zeros(1,length(signal));
    high_indexes = zeros(1,length(signal));
    for i = 1:length(signal)
        indexes = find(f(:,i) > 0.03 & f(:,i) < 0.5);
        low_indexes(i) = indexes(1);
        high_indexes(i) = indexes(end);
    end
    f_plot = f(min(low_indexes):max(high_indexes),:);
    Pb_plot = Pb(min(low_indexes):max(high_indexes),:);
    figure
    mesh(t_plot,f_plot,Pb_plot,'FaceColor','interp')
    xline(indexes_end_phase,'w--','LineWidth',0.8)
    xlabel('Time (samples)')
    ylabel('Frequency (Hz)')
    zlabel('PSD')
    title(sprintf('Window length: %d',winlen))
    ylim([0.04 0.4])
    grid on
    view(2)
    if saving_flag
        output_file = fullfile(figure_output_folder, ...
            strcat(output_label,'_param_',num2str(winlen),'.png'));
        saveas(gcf,output_file);
    end
end
end
