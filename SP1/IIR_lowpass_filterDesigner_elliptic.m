function [] = IIR_lowpass_filterDesigner_elliptic(y, Fs, MAT_file_name, play_filtered_output)
elliptic = load(MAT_file_name);
filtered_elliptic = filtfilt(elliptic.SOS, elliptic.G, y);

figure
spektrogram(filtered_elliptic, Fs);
title('Spectrogram of the filtered signal (IIR low-pass - elliptic)');

if play_filtered_output == true
    player = audioplayer(filtered_output,Fs);
    playblocking(player)
end
end

