%% 1. semestralni prace z predmetu ZSI
% Jan Burian

%%
clc
close all
clear all

%% Loading signal from .wav file
filename = "veta.wav";
[y, Fs] = audioread(filename);  
t = length(y) / Fs; % time of signal
%Ts = 1 / Fs; % sampling period
x = linspace(0, t, length(y));

figure 
plot(x, y);
title("Loaded signal");
xlabel("Time [s]");
ylabel("Amplitude");

% Spectrogram
figure
window = 1024;
spectrogram(y, window, [], [], Fs, 'yaxis');
title("Spectrogram of loaded signal");