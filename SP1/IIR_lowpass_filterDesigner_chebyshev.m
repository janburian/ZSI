function [] = IIR_lowpass_filterDesigner_chebyshev(y, Fs, MAT_file_name, play_filtered_output)
chebyshev = load(MAT_file_name);
filtered_chebyshev = filtfilt(chebyshev.SOS, chebyshev.G, y);

spektrogram(filtered_chebyshev, Fs);
title('Spectrogram of the filtered signal (IIR lowpass - chebyshev)');

if play_filtered_output == true
    player = audioplayer(filtered_output,Fs);
    playblocking(player)
end
end

