function [segments] = divide_signal_into_segments(signal, segment_length)

% Calculate the total number of segments
num_segments = floor(length(signal) / segment_length);

% Initialize the cell array to store the segments
segments = cell(num_segments, 1);

% Divide the signal into segments
for i = 1:num_segments
    start_idx = (i-1) * segment_length + 1;
    end_idx = min(start_idx + segment_length - 1, length(signal));
    segments{i} = signal(start_idx:end_idx);
end

end

