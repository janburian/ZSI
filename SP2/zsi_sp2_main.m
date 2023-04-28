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
title("Signal");
xlabel("Time [s]");
ylabel("Amplitude");





