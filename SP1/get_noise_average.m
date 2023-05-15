function [noise_average] = get_noise_average(fft_samples)
noise_sum = 0;

for j = 1:num_segments_noise
   noise_sum = noise_sum + abs(fft_samples{j}); 
end

noise_average = (1/num_segments_noise) * noise_sum;

end

