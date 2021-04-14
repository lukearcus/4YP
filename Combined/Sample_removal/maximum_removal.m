function [deltas_out,reg_y_out,removed] = maximum_removal(coord,deltas,reg_y)

    deltas_out = deltas;
    reg_y_out = reg_y;
    
    [~,y_ind] = max(coord.y_max);
    deltas_out(y_ind) = [];
    coord.num_deltas = size(deltas_out,2);
    coord.y_max(y_ind) = [];
    reg_y_out(y_ind) = [];
    removed =1;
    
    
end