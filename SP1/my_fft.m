function X = my_fft(x)
% Compute the FFT of a signal x using the Cooley-Tukey algorithm

% Check if the input signal length is a power of 2
N = length(x);
if log2(N) ~= floor(log2(N))
    error('Input signal length must be a power of 2');
end

% Compute the FFT recursively
if N == 1
    X = x;
else
    X_even = my_fft(x(1:2:N));
    X_odd = my_fft(x(2:2:N));
    twiddle = exp(-2j*pi*(0:N-1)/N);
    X = [X_even + twiddle(1:N/2).*X_odd, X_even + twiddle(N/2+1:N).*X_odd];
end
