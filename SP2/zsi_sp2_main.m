%% 2. semestralni prace z predmetu ZSI
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

%% Empirical mode decomposition
num_iterations = 12;
num_imf = 12;
S_criterion = 6;

signals_cell = cell(num_iterations,1);
IMF_cell = cell(num_iterations,1);

signals_cell{1} = y;

for i = 1:num_imf
    IMF_cell{i} = get_imf(signals_cell{i}, num_iterations, S_criterion);
    signals_cell{i+1} = signals_cell{i} - IMF_cell{i};
end

% Plot of IMFs
figure
for u = [1, 2, 4, 8, 12]
    hold on
    plot(x, IMF_cell{u});
end
ylim([-0.8 0.8])
xlabel('Time [s]')
ylabel('Amplitude [dB]')
title('IMF comparison')
legend('IMF 1','IMF 2','IMF 4','IMF 8','IMF 12')

% Plot of IMFs (zoomed)
figure
for u = [1, 2, 4, 8, 12]
    hold on
    plot(x, IMF_cell{u});
end
xlim([4 5.5])
ylim([-0.12 0.12])
xlabel('Time [s]')
ylabel('Amplitude [dB]')
title('IMF comparison (from 4 to 5.5 s)')
legend('IMF 1','IMF 2','IMF 4','IMF 8','IMF 12')

% Subplot of components
figure
sgtitle('Signal components')
for i = 1:num_imf
    subplot(12,1,i);
    plot(x, IMF_cell{i});
    xlim([4.2 4.6])
    %ylim([-0.2 0.2])
    ylabel(sprintf('c%d', i))
    set(gca,'XTick',[], 'YTick', [])
end

%% Amplitude spectrum
figure; 
hold on;
n = length(y); 
for j = [1, 2, 4, 8, 12]
    fft_out = fft(signals_cell{j}); 
    P2 = abs(fft_out) / n; 
    P1 = P2(1:n/2+1); 
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(n/2))/n;
    plot(f, P1) 
end

xlim([0 6500])
ylim([0 0.006])
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
title('Amplitude spectrum')
legend('IMF 1','IMF 2','IMF 4','IMF 8','IMF 12')

%% Instantaneous frequency
num_iterations = 1; 
num_used_imf = 3;
res_freq = cell(num_used_imf, 1);
S_criterion = 6; 
array_iterations = [10, 25, 100];

for i = array_iterations
    used_imf = cell(num_used_imf, 1);
    
    signals = cell(num_used_imf, 1);
    signals{1} = y; 
    
    for j = 1:num_used_imf
        used_imf{j} = get_imf(signals{j}, i, S_criterion);
        signals{j+1} = signals{j} - used_imf{j};
        hilb_imf = hilbert(used_imf{j}); % compute the Hilbert transform of the input IMF
        inst_phase = unwrap(angle(hilb_imf)); % instantaneous phase of the IMF as the angle of the Hilbert transform
        inst_freq = diff(inst_phase) / (2*pi/Fs); % instantaneous frequency of the IMF as the derivative of the instantaneous phase
        %inst_freq = Fs / (2*pi) * diff(inst_phase);
        res_freq{num_iterations, j} = inst_freq;
    end 
    num_iterations = num_iterations + 1; 
end

t = length(res_freq{1,1}) / Fs; % time of signal
x = linspace(0, t, length(res_freq{1,1})); % generating x vector

for i = 1:num_used_imf
    figure
    for j =1:num_used_imf
        subplot(3,1,j)
        plot(x, res_freq{i,j});
        ylabel('Frequency [Hz]')
        xlabel('Time [s]')
        title(sprintf('Instantaneous frequency IMF %d (%d iterations)', i, array_iterations(j)));
    end
end






