%% 1. semestralni prace z predmetu ZSI
% Jan Burian

%%
clc
close all
clear all

%% Loading signal from .wav file
filename = "veta.wav";
[y, Fs] = audioread(filename);  

num_samples = length(y); 
t = length(y) / Fs; % time of signal
%Ts = 1 / Fs; % sampling period

% Spectrogram
figure; 
time_domain_spectrogram(y, Fs); % parasitic frequency: 6000 Hz

%% FIR and IIR filters
%% FIR lowpass filter based on filterDesigner toolbox
% Application of FIR lowpass filter on signal
Num = load('FIR_lowpass_designer_toolbox.mat');
filtered_output = filter(Num.Num,1,y);
n = filtord(Num.Num, 1); % filter order

% Spectrogram
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; FIR lowpass filter (n = %d)', n))

% Play filtered signal
player = audioplayer(filtered_output,Fs);
playblocking(player)

%% FIR highpass filter based on filterDesigner toolbox
Num = load('FIR_highpass_designer_toolbox.mat');
filtered_output = filter(Num.Num,1,y);
n = filtord(Num.Num, 1); % filter order

% Application of FIR lowpass filter on signal
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; n = %d', n))

% Play filtered signal
player = audioplayer(filtered_output,Fs);
playblocking(player)

%% FIR lowpass filter
n    = 200;         % filter order
Fc   = 0.25;        % cutoff frequency (normalized)
flag = 'noscale';  % no normalization

% Windows
%win = rectwin(n);
%win = hann(n);
win = hamming(n); % best
%win = blackman(n);

% FIR_coeffs
b  = fir1(n-1, Fc, 'low', win, flag);    
fvtool(b,1)

% Filter application
filtered_output = filter(b,1,y);  

% Spectrogram
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; n = %d', n))

% Play filtered signal
player = audioplayer(filtered_output,Fs);
playblocking(player)

%% FIR band-pass filter
n  = 255;       % filter order for band-pass = n*2
Fc1  = 5700;    % 1st cutoff frequency [Hz]    
Fc2  = 6250;    % 2nd cutoff frequency [Hz]
flag = 'noscale';  

%win = rectwin(n);
%win = hann(n);
win = hamming(n);
%win = blackman(n);

b  = fir1(n-1, [Fc1 Fc2]/(Fs/2), 'stop', win, flag);  % getting coefs
%fvtool(b,1)

filtered_output = filter(b,1,y);    % application of the filter

% Spectrogram
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; n = %d', n))

% Play filtered signal
player = audioplayer(filtered_output,Fs);
playblocking(player)

%% IIR lowpass filter based on filterDesigner toolbox
% Chebyshev
chebyshev = load('IIR_lowpass_chebyshev_I.mat');
filtered_chebyshev = filtfilt(chebyshev.SOS, chebyshev.G, y);

spektrogram(filtered_chebyshev, Fs);
title('Spektrogram signalu po aplikaci IIR lowpass filtru');

player = audioplayer(filtered_output,Fs);
playblocking(player)

%% IIR band-stop filter based on filterDesigner toolbox
chebyshev = load('IIR_bandstop_chebyshev_II.mat');
filtered_chebyshev = filtfilt(chebyshev.SOS, chebyshev.G, y);

spektrogram(filtered_chebyshev, Fs);
title('Spektrogram signalu po aplikaci IIR band-stop filtru');

player = audioplayer(filtered_output,Fs);
playblocking(player)

%% IIR band-pass filter
n  = 1;       % filter order for band-pass = n*2
Fc1  = 5900;    % 1st cutoff frequency [Hz]    
Fc2  = 6100;    % 2nd cutoff frequency [Hz]

[B,A] = butter(n,[(Fc1/(0.5*Fs)) (Fc2/(0.5*Fs))],'stop');

filtered_output = filter(B,A,y);

% Spectrogram
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; n = %d', n))

% Play filtered signal
player = audioplayer(filtered_output,Fs);
playblocking(player)

%% Counting own filter (band-stop FIR filter)
% Define the filter specifications
Fc1 = 5400; % Lower cutoff frequency
Fc2 = 6500; % Upper cutoff frequency
N = 200; % Filter order

% Compute the filter coefficients
wc1 = 2*pi*Fc1/Fs; % Normalized lower stopband frequency
wc2 = 2*pi*Fc2/Fs; % Normalized upper stopband frequency
b = zeros(1, N+1); % Initialize the filter coefficients
for n = 0:N
    if n == N/2
        b(n+1) = 1 - (wc2-wc1)/pi;
    else
        b(n+1) = (sin(wc1*(n-N/2)) - sin(wc2*(n-N/2)))/(pi*(n-N/2));
    end
end

% Apply the filter to a signal
filtered_output = filter(b, 1, y); % Apply the filter to the signal

% Spectrogram
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; n = %d', n))

% Play filtered signal
player = audioplayer(filtered_output,Fs);
playblocking(player)

%%  



