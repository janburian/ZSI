clc
close all
clear all

%% Loading signal from .wav file
filename = "veta2.wav";
[y, Fs] = audioread(filename);  
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

max_indices = find(d_1 .* d_2 < 0 & d_1 > 0) + 1; % +1 jinak nevychazi
min_indices = find(d_1 .* d_2 < 0 & d_1 < 0) + 1; % +1 jinak nevychazi

extrema_max = zeros(length(max_indices), 2);
extrema_min = zeros(length(min_indices), 2);

for i=1:1:length(max_indices)
    idx_max = max_indices(i); 
    idx_min = min_indices(i); 
    
    extrema_max(i, 1) = y(idx_max);
    extrema_max(i, 2) = max_indices(i);
    
    extrema_min(i, 1) = y(idx_min);
    extrema_min(i, 2) = min_indices(i);
end

%% Zero crossings
num_zero_crossings = 0;
for j=1:1:length(y)-1
    if y(j) > 0 && y(j+1) < 0 || y(j) < 0 &&  y(j+1) > 0
        num_zero_crossings = num_zero_crossings + 1;
    end
end


%% Empirical mode decomposition
[imf,residual] = emd(y);




