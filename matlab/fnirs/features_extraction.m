function [features,diagnostics] = features_extraction( ...
        signal,time,i_phases,fs,gr)
%FEATURES_EXTRACTION Extract phase-level fNIRS concentration features.
%
% Inputs
% ------
% SIGNAL   : Hemoglobin concentration-change signal.
% TIME     : Time vector in seconds, same length as SIGNAL.
% I_PHASES : 2-by-N matrix of one-based inclusive phase limits.
% FS       : Sampling frequency in Hz.
% GR       : Logical flag. If true, diagnostic figures are generated.

% Feature names returned in the first column.
feature_names = ["min"; "max"; "mean"; "variance"; "std"; "skewness"; "kurtosis"; "power_band"; "max_frequency"; "median_frequency"; "spectral_entropy"; "mean_diff"; "peak_to_peak"; "slope"; "zero_crossing"; "polarity"; "entropy"; "auc"; "total_energy"];

% Set number of phases.
[~,number_phase] = size(i_phases);
 
% Keep names and numeric doubles in separate cell columns so feature values
% retain full precision.
features = cell(length(feature_names),number_phase+1);
features(:,1) = cellstr(feature_names);

% Apply a high-pass filter to remove slow concentration drift.
hp = designfilt( ...
    'highpassiir','PassbandFrequency',0.01, ...
    'StopbandFrequency',0.006,'PassbandRipple',0.5, ...
    'StopbandAttenuation',80,'SampleRate',fs, ...
    'DesignMethod','butter');

% Display the filter response.
if gr
    fvtool(hp);
end
% Signal filtered with a high-pass filter.
pad_seconds = 180;
pad_samples = round(pad_seconds * fs);
pad_samples = min(pad_samples, length(signal) - 2);

signal_filt = fnirs_zero_phase_sos( ...
    hp.Coefficients,signal,pad_samples);

diagnostics = struct( ...
    'FilteredSignal',signal_filt, ...
    'Time',time, ...
    'SamplingFrequency',fs, ...
    'PhaseIndexes',i_phases);

% Plot the signal before and after high-pass filtering.
if gr
    figure()
    plot(time, signal,'LineWidth',1.1)
    hold on
    plot(time, signal_filt,'LineWidth',1.1)
    xline(time(i_phases(1,:)),'--')
    xline(time(i_phases(2,:)),'--')
    axis tight
    hold off
    ylabel('Concentration change [mM]')
    xlabel('Time [s]')
    title('fNIRS signal - Before and after HPF')
    lgd = legend('Signal before HPF','Signal after HPF');
    lgd.FontSize = 12;
    lgd.FontWeight = 'bold';
end

% Iterate through the requested recording phases.
for i=1:size(i_phases,2)
    start_phases = i_phases(1,i);
    end_phases = i_phases(2,i);

    signal_phase = real(signal(start_phases:end_phases,1));
    signal_phase_filt = real(signal_filt(start_phases:end_phases,1));
    time_phase = time(start_phases:end_phases,1);

    % Calculate amplitude, time-domain, and spectral features.
    %% MEAN VALUE
    mean_feature = mean(signal_phase);
    
    %% ENTROPY
    if max(signal_phase) == min(signal_phase)
    entropy_feature = 0;
    else
        edges = linspace(min(signal_phase), max(signal_phase), 11);
        probability = histcounts(signal_phase, edges, 'Normalization', 'probability');
        probability = probability(probability > 0);
        entropy_feature = -sum(probability .* log2(probability));
    end

    %% AREA UNDER THE CURVE
    % Calculate area under the curve for the current phase.
    AUC = sum(0.5 * diff(time_phase) .* abs(signal_phase(1:end-1,1) + signal_phase(2:end,1)));
    AUC_feature = AUC/(time_phase(end,1)-time_phase(1,1));

    %% TOTAL ENERGY
    % Calculate total energy normalized by phase duration.
    total_energy_feature = sum(signal_phase.^2/(time_phase(end,1)-time_phase(1,1)));

    %% SLOPE
    P=polyfit(time_phase,signal_phase,1);
    % Extract the slope from the polynomial fit.
    slope_feature = P(1);
    y_fit = polyval(P,time_phase);

    if gr==1
        figure;
        plot(time_phase, signal_phase,'r', 'DisplayName', 'original signal');
        hold on;
        plot(time_phase, y_fit, '-','LineWidth',2,'DisplayName', 'Line of best fit');
        xlim([min(time_phase) max(time_phase)])
        slope_label = sprintf('Slope: %.2e', slope_feature);
        text('Position', [0.02, 0.9], 'String', slope_label, 'FontSize', 12, 'Units', 'normalized');
        hold off;
        xlabel('Time [s]')
        ylabel('Concentration change [mM]')
        legend('show');
    end

    %% Features calculated after detrending
    %% Statistical features
    %% MAX VALUE
    max_feature = max(signal_phase_filt);
    %% MIN VALUE
    min_feature = min(signal_phase_filt);
    %% KURTOSIS
    kurtosis_feature = kurtosis(signal_phase_filt,1);
    %% Skewness
    skewness_feature = skewness(signal_phase_filt,1);
    %% VARIANCE
    var_feature = var(signal_phase_filt);
    %% STANDARD DEVIATION
    std_feature = std(signal_phase_filt);

    %% TEMPORAL FEATURES
    %% ZERO CROSSING RATE
    % Count the fraction of samples crossing zero within the phase.
    N = size(signal_phase_filt,1);
    sample_sign = zeros(N,1);
    zcr = zeros(N,1);           % Vector to store zero crossing count result
    % Determine sample signs.
    for m=1:N
        if signal_phase_filt(m,1)>0
            sample_sign(m,1) = 1;
        elseif (signal_phase_filt(m,1)<0)
            sample_sign(m,1) = 0;
        end
    end
    % Loop for counting zero crossings
    for m = 1:N-1
        % Compare the signs of adjacent samples.
        delta_sign = sample_sign(m+1,1) - sample_sign(m,1);
        if delta_sign> 0 || delta_sign<0
            zcr(m,1) =1;   % Sign change, detect a zero crossing
        else
            zcr(m,1) =0;   % No sign change
        end
    end
    % Calculation of the average Zero Crossing Rate
    zero_crossing_feature = sum(zcr)/N;

    %% MEAN DIFFERENCE
    mean_difference_feature = mean(diff(signal_phase_filt));
    %% POLARITY
    % Polarity is the absolute ratio between the maximum and minimum values.
    max_value = max(signal_phase_filt);
    min_value = min(signal_phase_filt);
    polarity_feature = abs(max_value/min_value);


    %% SPECTRAL FEATURES
    %% FFT calculation
    % Remove the mean value from the windowed signal.
    signal_phase_no_mean = signal_phase_filt -mean(signal_phase_filt);
    % Calculate the length of the current phase window.
    N = length(signal_phase_no_mean);
    n = 2^nextpow2(N); % next power of 2 from the original signal length
    f_entropy = fs * (0:(n/2)) / n;
    L = (n / 2) + 1;

    % Calculate the FFT magnitude.
    fmag = abs(fft(signal_phase_no_mean,n));

    %% Spectral entropy
    % Calculate spectral entropy from the normalized power spectrum.
    power = fmag(1:L,1) .^ 2;

    if sum(power) == 0
        spectral_entropy_feature = 0;
    else
        % Normalize the power spectrum to obtain a probability distribution.
        prob = power / sum(power);
        prob = prob(prob ~= 0);

        % If probability is concentrated in one bin, entropy is zero.
        if numel(prob) == 1
            spectral_entropy_feature = 0;
        else
            spectral_entropy_feature = -sum(prob .* log2(prob)) / log2(numel(prob));
        end
    end

    %% Median frequency
    % Frequency where cumulative Fourier magnitude reaches 50%.
    cum_f_phase = cumsum(fmag(1:L,:),1);

    % Find the first index corresponding to 50% cumulative magnitude.
    ind_median = find(cum_f_phase > (cum_f_phase(end) * 0.50), 1);

    % Fallback to MATLAB's median index if no cumulative index is found.
    if isempty(ind_median)
        [~, ind_median] = median(cum_f_phase);
    end

    median_frequency_feature = f_entropy(ind_median);

    %% Welch PSD features
    % Estimate PSD over overlapping windows after removing each window mean.
    N=512/2;               % window points
    overlap = N/2;       % 50% window overlap
    % Calculate number of overlapping windows.
    window_flag=0;
    if mod(length(signal_phase_filt) ,N)~=0
        num_windows = fix(length(signal_phase_filt) / (overlap));
        window_flag=1;
    else
        num_windows = (length(signal_phase_filt) / (overlap))-1;
    end
    % Calculate PSD with pwelch for each window and average the estimates.
    Ptot = [];
    dim=[];
    f_tot=[];
    for z = 1 : num_windows
        overlap = N/2;
        start_idx = (z - 1) * (N-overlap) + 1;
        if  (window_flag == 1) && (z==num_windows)
            end_idx = length(signal_phase_filt);
        else
            end_idx = start_idx + N - 1;
        end

        window_phase= signal_phase_filt(start_idx:end_idx,1)-mean(signal_phase_filt(start_idx:end_idx,1),1);
        if length(window_phase)<N
            zero_new=N-length(window_phase);
            % Zero-padding.
            window_phase(length(window_phase)+1:length(window_phase)+zero_new,1) =zeros(zero_new,1);
        end
        dim =[dim; length(window_phase)];

        overlap=0;
        nfft = 2048;
        window = hamming(dim(z,1));
        [pxx, f] = pwelch(window_phase, window, overlap, 2*nfft, fs);

        Ptot = [Ptot, pxx];
        f_tot = [f_tot, f];
    end
    power=mean(Ptot(:,1:end),2);
    freq = f_tot(:,1);

    if gr==1
        figure
        plot(freq, power/max(power), 'k', 'LineWidth', 1.2);
        xlabel('Frequency (Hz)');
        ylabel('Power Spectral Density');
        title(['PSD SIGNAL PHASE:',num2str(i)]);
        grid on;
        xlim([0 0.8])
    end

    %% Max frequency
    % Frequency corresponding to the maximum PSD peak.
    [~,index_max_power] = max(power);
    max_frequency_feature = freq(index_max_power);

    %% Power bandwidth
    % Width of the frequency band containing 95% of the PSD power.

    if sum(power) == 0
        power_bandwidth_feature = 0;
    else

        % Compute lower and upper power-bandwidth limits.
        cum_power = cumsum(power);
        f_lower = freq(find(cum_power >= (cum_power(end) * 0.95), 1));

        % Calculate inverse cumulative power.
        cum_power_inv = cumsum(flipud(power));
        index_frequency = abs(find(cum_power_inv >= (cum_power(end) * 0.95), 1) - size(power,1) + 1);
        if index_frequency == 0
            index_frequency = 1;
        end
        f_upper = freq(index_frequency);

        power_bandwidth_feature = abs(f_upper - f_lower);
    end

    % Smooth the phase signal before estimating peak-to-peak distance.
    cutoff_frequency = 0.08;
    filter_order = 5;

    % Butterworth low-pass filter.
    [b, a] = butter(filter_order, cutoff_frequency / (fs/2), 'low');
    [sos,scale] = tf2sos(b,a);
    sos(1,1:3) = scale*sos(1,1:3);

    % Display filter response.
    if gr==1
        freqz(b, a, 1024, fs);
    end
    % Smoothed phase signal.
    smoothing_signal_phase = fnirs_zero_phase_sos( ...
        sos,signal_phase_filt,3*filter_order);

    %% Peak-to-peak distance
    % Calculate the maximum distance between a local maximum and the next
    % local minimum.

    % Calculate the first derivative.
    d_signal = diff(smoothing_signal_phase);

    % Locate local maxima and minima in the smoothed phase signal.
    max_local_index = find(d_signal(1:end-1,1)>0 & d_signal(2:end,1)<=0)+1;
    min_local_index = find(d_signal(1:end-1,1)<0 & d_signal(2:end,1)>=0)+1;

    % If either local-extrema vector is empty, set the distance to zero.
    if isempty(max_local_index) || isempty(min_local_index)
        peaktopeakdistance_feature = 0;
    else
        % For a single maximum/minimum pair, require the minimum to follow
        % the maximum.
        if length(min_local_index)==length(max_local_index) && length(min_local_index)==1
            if min_local_index(1,1)< max_local_index(1,1)
                peaktopeakdistance_feature = 0;
            else
                peaktopeakdistance_feature = abs(smoothing_signal_phase(max_local_index(1,1),1)-smoothing_signal_phase(min_local_index(1,1),1));
            end
        else
            % Remove a leading minimum that has no preceding local maximum.
            if min_local_index(1,1)< max_local_index(1,1)
                min_local_index = min_local_index(2:end,1);
            end
            % Remove a trailing maximum that has no following local minimum.
            if max_local_index(end,1)> min_local_index(end,1)
                max_local_index = max_local_index(1:end-1,1);
            end
            if isempty(max_local_index) || isempty(min_local_index)
                peaktopeakdistance_feature = 0;
            else
                % Calculate distance from each local maximum to the next
                % local minimum.
                peaktopeakdistance=zeros(length(max_local_index),1);
                for k=1:length(max_local_index)
                    peaktopeakdistance(k,1) = abs(smoothing_signal_phase(max_local_index(k,1),1)-smoothing_signal_phase(min_local_index(k,1),1));
                end
                % Select the maximum distance.
                peaktopeakdistance_feature = max(peaktopeakdistance);
            end
        end
    end
    % Plot peak-to-peak feature diagnostics.
    if gr
        figure()
        plot(time_phase,signal_phase_filt)
        hold on
        plot(time_phase,smoothing_signal_phase,'LineWidth',1.2);
        plot(time_phase(max_local_index,1),smoothing_signal_phase(max_local_index,1),'ro','LineWidth',2)
        plot(time_phase(min_local_index,1),smoothing_signal_phase(min_local_index,1),'bo','LineWidth',2)
        hold off
        grid on
        xlim([min(time_phase) max(time_phase)])
        legend('original signal','filtered signal with LPF','max local','min local')
        xlabel('Time [s]')
        ylabel('Concentration change [mM]')
    end
    phase_feature_values = [min_feature; max_feature; mean_feature; ...
        var_feature; std_feature; skewness_feature; kurtosis_feature; ...
        power_bandwidth_feature; max_frequency_feature; ...
        median_frequency_feature; spectral_entropy_feature; ...
        mean_difference_feature; peaktopeakdistance_feature; ...
        slope_feature; zero_crossing_feature; polarity_feature; ...
        entropy_feature; AUC_feature; total_energy_feature];
    features(:,i+1) = num2cell(phase_feature_values);
end
