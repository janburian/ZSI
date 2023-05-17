function [X] = my_fft(x)
% Compute the FFT of a signal x using Cooley–Tukey FFT algorithm

% Check if the input signal length is a power of 2
N = length(x);
if log2(N) ~= floor(log2(N))
    error('Input signal length must be a power of 2.');
end

% If signal has only one element
if N == 1
    X = x;
else
    % Extracting even and odd elements
    even_elements = x(1:2:N);
    odd_elements = x(2:2:N);
    
    % Compute the FFT recursively 
    X_even = my_fft(even_elements);
    X_odd = my_fft(odd_elements);
    
    % Complex part
    twiddle_factor = exp(-2j*pi*(0:N-1)/N);
    X = [X_even + twiddle_factor(1:N/2).*X_odd, X_even + twiddle_factor(N/2+1:N).*X_odd];
end
