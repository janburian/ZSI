function [y_resampled] = resample_signal(y, fs_orig, fs_new)
% Resample a signal y from fs_orig to fs_new using linear interpolation

% Compute the resampling factor
r = fs_new / fs_orig;

% Compute the resampled signal length
N_resampled = round(length(y) * r);

% Initialize the resampled signal
y_resampled = zeros(1, N_resampled);

% Resample the signal using linear interpolation
for n_resampled = 1:N_resampled
    n_orig = (n_resampled - 0.5) / r + 0.5;
    n1 = floor(n_orig);
    n2 = ceil(n_orig);
    if n2 > length(y)
        y_resampled(n_resampled) = y(end);
    else
        y_resampled(n_resampled) = y(n1) + (n_orig - n1) * (y(n2) - y(n1));
    end
end

end