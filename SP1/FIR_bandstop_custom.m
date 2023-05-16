function [] = FIR_bandstop_custom(y, Fs, N, Fc1, Fc2, play_filtered_output)
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
title(sprintf('Spectrogram of the filtered signal; n = %d', N))

% Play filtered signal
if play_filtered_output == true
    player = audioplayer(filtered_output,Fs);
    playblocking(player)
end
end

