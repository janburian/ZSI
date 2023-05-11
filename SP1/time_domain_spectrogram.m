function time_domain_spectrogram(y,fs)
% Parameters
window = round(0.03*fs); 
overlap = round(0.5*window);
FFT_resolution = window;

% Spectrogram
[y1,f1,t1,p1] = spectrogram(y,window,overlap,FFT_resolution,fs);

% Plotting
subplot(2,1,1);
plot((1:1:length(y))/fs,y);
axis tight; 
xlabel('Time [s]');
ylabel('Amplitude');
title('Loaded original signal in time domain');

subplot(2,1,2);
surf(t1,f1,log10(abs(p1)),'EdgeColor','none');
axis xy; 
axis tight; 
colormap(jet); 
view(0,90);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('Spectrogram of the original signal');
ylim([0 15000])



