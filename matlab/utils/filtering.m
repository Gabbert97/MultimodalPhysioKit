function sig_f = filtering(signal,fs,ft_high,ft_low,ord_high,ord_low,figures_1_0)
%FILTERING Apply zero-phase high-pass, low-pass, or band-pass filtering.
%
% Inputs
% ------
% SIGNAL      : Input signal vector.
% FS          : Sampling frequency in Hz.
% FT_HIGH     : High-pass cutoff frequency in Hz. Use 0 to disable.
% FT_LOW      : Low-pass cutoff frequency in Hz. Use 0 to disable.
% ORD_HIGH    : High-pass Butterworth filter order.
% ORD_LOW     : Low-pass Butterworth filter order.
% FIGURES_1_0 : If 1, display filter responses.
%
% Output
% ------
% SIG_F       : Filtered signal.

% Nyquist frequency
FNy = fs/2;
% HPF
Wn_high = ft_high/FNy;
% LPF
Wn_low = ft_low/FNy;

if ft_high > 0 && ft_low > 0
    %% HighPass Filter
    [b,a] = butter(ord_high,Wn_high,'high');
    if figures_1_0 == 1
        figure('Name','FPA')
        freqz(b,a,1024,fs);
    end
    sig_fpa = filtfilt(b,a,double(signal));
    %% LowPass Filter
    [b,a] = butter(ord_low,Wn_low,'low');
    if figures_1_0 == 1
        figure('Name','FPB')
        freqz(b,a,1024,fs);
    end
    sig_f = filtfilt(b,a,sig_fpa);
elseif ft_low == 0
    %% HighPass Filter
    [b,a] = butter(ord_high,Wn_high,'high');
    if figures_1_0 == 1
        figure('Name','FPA')
        freqz(b,a,1024,fs);
    end
    sig_f = filtfilt(b,a,double(signal));
elseif ft_high == 0
    %% LowPass Filter
    [b,a] = butter(ord_low,Wn_low,'low');
    if figures_1_0 == 1
       	figure('Name','FPB')
        freqz(b,a,1024,fs);
    end
    sig_f = filtfilt(b,a,double(signal));
end
end
