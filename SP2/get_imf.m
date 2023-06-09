% IMF
function [signal_out] = get_imf(signal, num_iterations, S_criterion)
    mean_value = 0;
    num_min_values = 0; 
    num_max_values = 0; 
    S_criterion_temp = 0; 
    num_zero_crossings_temp = 0;
    
    
    for i = 1:num_iterations
        [extremes_min, extremes_max] = get_local_extremes(signal);
        num_zero_crossings = get_zero_crossings(signal);
        
        if (num_min_values == length(extremes_min(:,1)) && num_max_values == length(extremes_max(:,1))...
            && num_zero_crossings_temp == num_zero_crossings)
            S_criterion_temp = S_criterion_temp + 1;
        else 
            S_criterion_temp = 0;
        end
        
        if (length(extremes_min(1,:)) < 1 || length(extremes_max(1,:)) < 1 || ...
                S_criterion_temp == S_criterion)
            break
        end
        
        % Interpolation
        min_interpolation = spline(extremes_min(:,1), extremes_min(:,2), 1:1:length(signal)); 
        max_interpolation = spline(extremes_max(:,1), extremes_max(:,2), 1:1:length(signal)); 
        
        num_min_values = length(extremes_min(1,:));
        num_max_values = length(extremes_max(1,:));
        num_zero_crossings_temp = num_zero_crossings; 
        
        % Mean value of envelope
        mean_value = (min_interpolation + max_interpolation) / 2; 
        
        % Component
        signal = signal - mean_value'; 
    end
    
    signal_out = signal; 
end

