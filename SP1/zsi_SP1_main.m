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
title(sprintf('Spectrogram of the filtered signal; n = %d', n))

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
n    = 150;         % filter order
Fc   = 0.25;        % cutoff frequency (normalized)
flag = 'noscale';  % no normalization

% Windows
%win = rectwin(n);
%win = hann(n);
win = hamming(n);
%win = blackman(n);

% FIR_coeffs
b  = fir1(n-1, Fc, 'low', win, flag);    
%fvtool(b,1)

% Filter application
filtered_output = filter(b,1,y);  

% Spectrogram
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; n = %d', n))

% Play filtered signal
player = audioplayer(filtered_output,Fs);
playblocking(player)

