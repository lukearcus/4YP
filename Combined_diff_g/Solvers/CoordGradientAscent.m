function [y_max,flag, cost] =CoordGradientAscent(coord,~,~,deltas,~,step_size,~)
    
    tot_demand = coord.tot_demand;
    grad = zeros(coord.num_deltas,1);
    for i = 1:coord.num_deltas
       grad(i) = tot_demand'*(diag(deltas(i).pr(:,1))*tot_demand + deltas(i).pr(:,2)); 
    end
    new_y = coord.y_max + step_size * grad;
    y_max = coord.poly.project_onto(new_y);
    cost = 0;
    for i = 1:coord.num_deltas
        cost = cost + new_y(i)*(tot_demand'*(diag(deltas(i).pr(:,1))*tot_demand + deltas(i).pr(:,2)));
    end
    flag = 1;
end