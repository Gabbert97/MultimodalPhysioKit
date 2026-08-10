function [S,t,f] = STFT(x,nfft,fs,winlen,window)
%STFT Compute a centered short-time Fourier transform.
%
% Inputs
% ------
% X      : Input signal vector.
% NFFT   : FFT length.
% FS     : Sampling frequency in Hz.
% WINLEN : Window length in samples.
% WINDOW : Window vector of length WINLEN.
%
% Outputs
% -------
% S : Complex STFT matrix for non-negative frequencies.
% T : Time vector in seconds.
% F : Frequency vector in Hz.

t=0:1/fs:((length(x)-1)/fs);
f=0:fs/nfft:(fs/2);

S=zeros(winlen,length(x));
if size(x,2)==1
    x=x.';
end

if rem(winlen,2) == 0
    x1=[zeros(1,winlen/2-1) x zeros(1,winlen/2)];
else
    x1=[zeros(1,(winlen-1)/2) x zeros(1,(winlen-1)/2)];   
end

for i=1:length(x)   
   x2 = x1(i:winlen-1+i) - mean(x1(i:winlen-1+i));
   S(:,i)=window.*x2';
end
    
S = fft(S,nfft,1);
S = S(1:(nfft/2)+1,:);
