function [] = FIR_lowpass_filterDesigner(y, Fs, MAT_file_name, play_filtered_output)
% Application of FIR lowpass filter on signal
Num = load(MAT_file_name);
filtered_output = filter(Num.Num,1,y);
n = filtord(Num.Num, 1); % filter order

fvtool(Num.Num, 1);

% Spectrogram
figure
spektrogram(filtered_output, Fs);
title(sprintf('Spectrogram of the filtered signal; FIR lowpass filter (n = %d)', n))

% Play filtered signal
if play_filtered_output == true
    player = audioplayer(filtered_output,Fs);
    playblocking(player)
end
end

