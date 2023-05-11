%% Identification of local extremes
function [extremes_min, extremes_max] = get_local_extremes(y)
    d = diff(y); 
    d_1 = d(1:end-1);
    d_2 = d(2:end);

    max_indices = find(d_1 .* d_2 < 0 & d_1 > 0) + 1; % +1 jinak nevychazi
    min_indices = find(d_1 .* d_2 < 0 & d_1 < 0) + 1; % +1 jinak nevychazi

    extremes_max = zeros(length(max_indices), 2);
    extremes_min = zeros(length(min_indices), 2);

    for i=1:1:length(max_indices)
        idx_max = max_indices(i); 

        extremes_max(i, 1) = max_indices(i);
        extremes_max(i, 2) = y(idx_max);
    end
    
    for i=1:1:length(min_indices)
        idx_min = min_indices(i); 

        extremes_min(i, 1) = min_indices(i);
        extremes_min(i, 2) = y(idx_min);
    end
end

