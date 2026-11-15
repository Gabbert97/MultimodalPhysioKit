function y = antialias_downsample(x,down,fsNew)
%ANTIALIAS_DOWNSAMPLE Low-pass filter, then keep every DOWN-th sample.
%
% Inputs
% ------
% X     : Signal vector sampled at DOWN*FSNEW Hz.
% DOWN  : Positive integer decimation factor.
% FSNEW : Target sampling frequency [Hz].
%
% Output
% ------
% Y     : X(1:DOWN:end) after a zero-phase 4th-order Butterworth low-pass
%         at 0.4*FSNEW, so content above the new Nyquist frequency cannot
%         alias into the analysis band. Orientation of X is preserved.
%         With DOWN == 1, X is returned unfiltered.

if down == 1
    y = double(x);
    return
end
wasRow = isrow(x);
fs = fsNew * down;
[z,p,k] = butter(4,(0.4*fsNew)/(fs/2),'low');
[sos,g] = zp2sos(z,p,k);
filtered = filtfilt(sos,g,double(x(:)));
y = filtered(1:down:end);
if wasRow
    y = y.';
end
end
