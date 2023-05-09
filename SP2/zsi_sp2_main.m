%% 2. semestralni prace z predmetu ZSI
% Jan Burian

%%
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
% d = diff(y); 
% d_1 = d(1:end-1);
% d_2 = d(2:end);
% 
% max_indices = find(d_1 .* d_2 < 0 & d_1 > 0) + 1; % +1 jinak nevychazi
% min_indices = find(d_1 .* d_2 < 0 & d_1 < 0) + 1; % +1 jinak nevychazi
% 
% extrema_max = zeros(length(max_indices), 2);
% extrema_min = zeros(length(min_indices), 2);
% 
% for i=1:1:length(max_indices)
%     idx_max = max_indices(i); 
%     idx_min = min_indices(i); 
%     
%     extrema_max(i, 1) = y(idx_max);
%     extrema_max(i, 2) = max_indices(i);
%     
%     extrema_min(i, 1) = y(idx_min);
%     extrema_min(i, 2) = min_indices(i);
% end

[extrema_min, extrema_max] = get_local_extrema(y);


%% Zero crossings
% num_zero_crossings = 0;
% for j=1:1:length(y)-1
%     if y(j) > 0 && y(j+1) < 0 || y(j) < 0 &&  y(j+1) > 0
%         num_zero_crossings = num_zero_crossings + 1;
%     end
% end

num_zero_crossings = get_zero_crossings(y); 

%% Empirical mode decomposition
num_iterations = 12;
S_criterion = 6;

signals_cell = cell(num_iterations,1);
IMF_cell = cell(num_iterations,1);

signals_cell{1} = y;

for i = 1:12
    IMF_cell{i} = get_imf(signals_cell{i}, num_iterations, S_criterion);
    signals_cell{i+1} = signals_cell{i} - IMF_cell{i};
end

figure
hold on
plot(x, IMF_cell{1});
plot(x, IMF_cell{2});
plot(x, IMF_cell{4});
plot(x, IMF_cell{6});
plot(x, IMF_cell{10});
ylim([-0.8 0.8])
xlabel('Time [s]')
ylabel('Amplitude [dB]')
title('IMF comparison')
legend('1.IMF','2.IMF','4.IMF','6.IMF','10.IMF')

%% Amplitude spectrum
index = 0;
figure; 
hold on;
n = length(y); 
for j = [1, 2, 4, 6, 10]
    fft_out = fft(signals_cell{j}); 
    P2 = abs(fft_out) / n; 
    P1 = P2(1:n/2+1); 
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(n/2))/n;
    plot(f, P1) 
end

xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
title('Amplitude spectrum')

%% 






