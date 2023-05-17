function [] = FIR_bandstop(y, Fs, n, Fc1, Fc2, flag, window, play_filtered_output)
b  = fir1(n-1, [Fc1 Fc2]/(Fs/2), 'stop', window, flag);  % Getting coefs
fvtool(b,1)

filtered_output = filter(b,1,y);    % Application of the filter

% Spectrogram
figure
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; FIR band-stop filter (n = %d)', n))

% Play filtered signal
if play_filtered_output == true
    player = audioplayer(filtered_output,Fs);
    playblocking(player)
end
end

