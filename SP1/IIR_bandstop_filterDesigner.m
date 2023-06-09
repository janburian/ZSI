function [] = IIR_bandstop_filterDesigner(y, Fs, MAT_file_name, play_filtered_output)
chebyshev_II = load(MAT_file_name);
filtered_chebyshev_II = filtfilt(chebyshev_II.SOS, chebyshev_II.G, y);

figure
spektrogram(filtered_chebyshev_II, Fs);
title('Spectrogram of the filtered signal; IIR band-stop filter');

if play_filtered_output == true
    player = audioplayer(filtered_output,Fs);
    playblocking(player)
end
end

