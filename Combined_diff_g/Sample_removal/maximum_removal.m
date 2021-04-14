function [deltas_out,reg_y_out,removed] = maximum_removal(coord,deltas,y_val,reg_y)

    deltas_out = deltas;
    reg_y_out = reg_y;
    
    [~,y_ind] = max(y_val);
    deltas_out(y_ind) = [];
    for i = 1:size(coord,2)
        coord(i).num_deltas = size(deltas_out,2);
        coord(i).y_max(y_ind) = [];
    end
    reg_y_out(:,y_ind) = [];
    removed =1;
    
    
end