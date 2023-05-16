%% 1. semestralni prace z predmetu ZSI
% Jan Burian

%%
clc
close all
clear all

%% Loading signal from .wav file
filename = 'veta.wav';
[y, Fs] = audioread(filename);  

num_samples = length(y); 
t = length(y) / Fs; % time of signal
%Ts = 1 / Fs; % sampling period

% Spectrogram
time_domain_spectrogram(y, Fs); % parasitic frequency: 6000 Hz
spektrogram(y, Fs);
title('Spectrogram of the original signal')

%% 1. FIR and IIR filters
%% FIR lowpass filter based on filterDesigner toolbox
FIR_lowpass_filterDesigner(y, Fs, './MAT_files/FIR_lowpass_designer_toolbox.mat', false)

%% FIR lowpass filter
n    = 300;         % filter order
Fc   = 0.25;        % cutoff frequency (normalized) % 5900 / (0.5 * Fs)
flag = 'noscale';  % no normalization

% Windows
%win = rectwin(n);
%win = hann(n);
window = hamming(n); % best
%win = blackman(n);

FIR_lowpass(y, Fs, n, Fc, flag, window, false); 

%% FIR band-stop filter
n  = 255;       % filter order for band-pass = n*2
Fc1  = 5600;    % 1st cutoff frequency [Hz]    
Fc2  = 6250;    % 2nd cutoff frequency [Hz]
flag = 'noscale';  

%win = rectwin(n);
%win = hann(n);
window = hamming(n);
%win = blackman(n);

FIR_bandstop(y, Fs, n, Fc1, Fc2, flag, window, false)

%% IIR lowpass filter based on filterDesigner toolbox
% Chebyshev
IIR_lowpass_filterDesigner_chebyshev(y, Fs, './MAT_files/IIR_lowpass_chebyshev_I.mat', false)

% Elliptic
IIR_lowpass_filterDesigner_elliptic(y, Fs, './MAT_files/IIR_lowpass_elliptic.mat', false)

%% IIR band-stop filter based on filterDesigner toolbox
IIR_bandstop_filterDesigner(y, Fs, './MAT_files/IIR_bandstop_chebyshev_II.mat', false)

%% IIR band-stop filter
n  = 1;       % filter order for band-pass = n*2
Fc1  = 5900;    % 1st cutoff frequency [Hz]    
Fc2  = 6100;    % 2nd cutoff frequency [Hz]

IIR_bandstop(y, Fs, n, Fc1, Fc2, false)

%% Custom FIR low pass filter
% Define the filter specifications
Fc = 5800; % Cut-off frequency
N = 1000; % Filter order

FIR_lowpass_custom(y, Fs, N, Fc, false)

%% Counting own filter (band-stop FIR filter)
% Define the filter specifications
Fc1 = 5700; % Lower cutoff frequency
Fc2 = 6100; % Upper cutoff frequency
N = 1500; % Filter order

FIR_bandstop_custom(y, Fs, N, Fc1, Fc2, false)

%% 2. Changing sampling frequency
target_Fs = 8000;

% Application of IIR low pass filter (antialising filter) to prevent aliasing effect       
% getting rid of the frequencies > 4000 

chebyshev_resampling = load('./MAT_files/IIR_resampling.mat');
filtered_chebyshev_resampling = filtfilt(chebyshev_resampling.SOS, chebyshev_resampling.G, y);

% Spectrogram
spektrogram(filtered_chebyshev_resampling, Fs);
title('Spectrogram of the filtered signal; IIR low pass filter')

resampled_signal = resample_signal(filtered_chebyshev_resampling, Fs, target_Fs);
spektrogram(resampled_signal, target_Fs);

output_filename = 'resampled_sentence.wav';
audiowrite(output_filename, resampled_signal, target_Fs);

%% 3. Removing additive noise
% Parameters
N = 8;
segment_length = 2^N; % m
alpha = 0.85;

% Dividing signal into segments 
segments = divide_signal_into_segments(resampled_signal, segment_length);
num_segments = length(segments);

% Computing the FFT of the resampled signal
fft_samples = cell(num_segments, 1);
for i = 1:num_segments
    fft_samples{i} = my_fft(segments{i});
end

% Parameters of the resampled signal
resampled_signal_length = length(resampled_signal);
resampled_signal_time = resampled_signal_length / target_Fs;

% Parameters of the noise which we want to remove
noise_time = 3;  % 3 seconds with noise and without speech
num_samples = noise_time * (resampled_signal_length / resampled_signal_time);
num_segments_noise = floor(num_samples / segment_length);
noise_average = get_noise_average(num_segments_noise, fft_samples);

% Substraction of the noise
fft_edited = cell(num_segments,1);
for i=1:num_segments
    substraction = abs(fft_samples{i}) - alpha*noise_average;
    if any(substraction < 0)
        substraction(substraction < 0) = 0;
    end
    fft_edited{i} = ifft(substraction .* exp(1i*angle(fft_samples{i}))); 
end

% Cell to vector
signal_without_noise_vec = real([fft_edited{:}]);

% Saving the signal without noise      
audiowrite('signal_without_noise.wav', signal_without_noise_vec, target_Fs);

spektrogram(signal_without_noise_vec, target_Fs);
title('Spectrogram of signal after noise removal');

subplot(2,1,1);
spektrogram(resampled_signal, target_Fs);
title('Spectrogram of the resampled signal')
subplot(2,1,2);
spektrogram(signal_without_noise_vec, target_Fs);
title('Spectrogram of signal after noise removal')







