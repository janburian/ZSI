function [] = FIR_lowpass(y, Fs, n, Fc, flag, window, play_filtered_output)
% FIR lowpass filter

% FIR_coeffs
b  = fir1(n-1, Fc, 'low', window, flag);    
fvtool(b,1)

% Filter application
filtered_output = filter(b,1,y);  

% Spectrogram
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; n = %d', n))

% Play filtered signal
if play_filtered_output == true
    player = audioplayer(filtered_output,Fs);
    playblocking(player)
end
end

