clc
close all
clear all
%% Loading signal from .wav file
filename = "veta.wav";
[y,Fs] = audioread(filename);  
t = length(y) / Fs; % time of signal
%Ts = 1 / Fs; % perioda vzorkovani
x = linspace(0, t, length(y));

figure; 
plot(x, y);
title("Loaded signal");
xlabel("Time [s]");
ylabel("Amplitude");

%% Identification of local extrema
d = diff(y); 
d_1 = d(1:end-1);
d_2 = d(2:end);

max_indices = find(d_1 .* d_2 < 0 & d_1 > 0);
min_indices = find(d_1 .* d_2 < 0 & d_1 < 0);

%% Empirical mode decomposition
[imf,residual] = emd(y);




