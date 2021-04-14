function [deltas_out,reg_y_out,removed] = all_active_removal(coord,deltas,reg_y,number_left)

    removed=0;
    threshold = 1e-5;
    deltas_out = deltas;
    reg_y_out = reg_y;
    
    while number_left > 0
        [y_val,y_ind] = max(coord.y_max);
        if y_val < threshold
           break 
        end
        deltas_out(y_ind) = [];
        coord.y_max(y_ind) = [];
        reg_y_out(y_ind) = [];
        number_left = number_left-1;
        removed = removed +1;
    end
    coord.num_deltas = size(deltas_out,2);
    

end