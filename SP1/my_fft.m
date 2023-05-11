function [X] = my_fft(x)
% Compute the FFT of a signal x

% Check if the input signal length is a power of 2
N = length(x);
if log2(N) ~= floor(log2(N))
    error('Input signal length must be a power of 2');
end

% Initialize the FFT coefficients
W = exp(-2j*pi*(0:N/2-1)/N);

% Perform the FFT
for m = 1:log2(N)
    for k = 1:2^(m-1)
        for n = 1:2^(log2(N)-m+1)
            idx1 = (k-1)*2^m + n;
            idx2 = idx1 - 2^(m-1);
            X(idx1) = x(idx1) + W(n*(2^(m-1))) * x(idx2);
            X(idx2) = x(idx1) - W(n*(2^(m-1))) * x(idx2);
        end
    end
    x = X; % Update the input signal for the next iteration
end
    
end
