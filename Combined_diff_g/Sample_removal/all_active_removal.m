function [deltas_out,reg_y_out,removed] = all_active_removal(coord,deltas,y_vals,reg_y,number_left)

    removed=0;
    threshold = 1e-5;
    deltas_out = deltas;
    reg_y_out = reg_y;
    
    while number_left > 0
        [y_val,y_ind] = max(y_vals);
        if y_val < threshold
           break 
        end
        deltas_out(y_ind) = [];
        y_vals(y_ind) = [];
        for i = 1:size(coord,2)
            coord(i).y_max(y_ind) = [];
            
        end
        reg_y_out(:,y_ind) = [];
        number_left = number_left-1;
        removed = removed +1;
    end
    for i = 1:size(coord,2)
        coord(i).num_deltas = size(deltas_out,2);
    end
    

end