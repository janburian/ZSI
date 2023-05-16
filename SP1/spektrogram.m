function spektrogram(y,fs)

% Plot signal spectrum

% Parameters
window = round(0.03*fs); 
overlap = round(0.5*window);
FFT_resolution = window;

% Spectrogram
[y1,f1,t1,p1] = spectrogram(y,window,overlap,FFT_resolution,fs);

surf(t1,f1,log10(abs(p1)),'EdgeColor','none');
axis xy; 
axis tight; 
colormap(jet); 
view(0,90);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
ylim([0 15000])


