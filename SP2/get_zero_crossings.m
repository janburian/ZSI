%% Zero crossings
function [num_zero_crossings] = get_zero_crossings(y)
    num_zero_crossings = 0;
    for j=1:1:length(y)-1
        if y(j) > 0 && y(j+1) < 0 || y(j) < 0 &&  y(j+1) > 0
            num_zero_crossings = num_zero_crossings + 1;
        end
    end
end

