function [] = my_FIR_lowpass(y, Fs, N, Fc, play_filtered_output)
% Compute the filter coefficients
wc = 2*pi*Fc/Fs; % Normalized cut-off frequency
b = zeros(1, N+1); % Initialize the filter coefficients
for n = 0:N
    if n == N/2
        b(n+1) = (1/pi) * wc;
    else
        b(n+1) = (1/(pi*(n-N/2))) * sin(wc*(n-N/2));
    end
end

% Apply the filter to a signal
filtered_output = filter(b, 1, y); % Apply the filter to the signal

% Spectrogram
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; n = %d', N))

% Play filtered signal
if play_filtered_output == true
    player = audioplayer(filtered_output,Fs);
    playblocking(player)
end
end

