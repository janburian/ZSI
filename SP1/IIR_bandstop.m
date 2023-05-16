function [] = IIR_bandstop(y, Fs, n, Fc1, Fc2, play_filtered_output)
[B,A] = butter(n,[(Fc1/(0.5*Fs)) (Fc2/(0.5*Fs))],'stop');

filtered_output = filter(B,A,y);

% Spectrogram
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; n = %d', n))

% Play filtered signal
if play_filtered_output == true
    player = audioplayer(filtered_output,Fs);
    playblocking(player)
end
end

