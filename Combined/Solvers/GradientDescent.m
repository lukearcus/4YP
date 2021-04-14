function [new_schedule,flag,f_cost,x_hat] = GradientDescent(ev,ev_fun_inputs)
    coord = ev_fun_inputs.coord;
    deltas = ev_fun_inputs.deltas;
    poly = ev.poly;
    step_size = ev_fun_inputs.coeff;
    t = ev_fun_inputs.t;
    
    tot_demand = coord.tot_demand;
    x_meno_i = tot_demand - ev.schedule;
    b_0 = zeros(size(tot_demand));
    curr_sched = ev.schedule;
    A_0 = diag(coord.A_0);
    
    f_gradient = 2*A_0*ev.schedule + A_0*x_meno_i+b_0;
    A_deltas = zeros(size(A_0));
    b_deltas = zeros(size(b_0));
     for i = 1:size(coord.y_max,2)
         A_deltas = A_deltas + coord.y_max(i) * diag(deltas(i).pr(:,1));
         b_deltas = b_deltas + coord.y_max(i) * (deltas(i).pr(:,2));
     end
     

    g_gradient = 2*A_deltas*tot_demand+b_deltas;
    
    gradient = f_gradient+g_gradient/coord.num_ag;
    sched = curr_sched - step_size*gradient;
    
    new_schedule = poly.project_onto(sched);
    flag = 1;
    
    %cost should be defined using x_hat
    x_hat = (ev.x_hat*(t-1) + new_schedule)/t;
    %x_hat = ((ev.x_hat*ev.eta_sum) + new_schedule*step_size)/(ev.eta_sum + step_size);
    %x_hat = poly.project_onto(x_hat); %to ensure x_hat remains in feasible set

    new_tot = coord.tot_demand_hat - ev.x_hat + x_hat;
    g_vals = zeros(1,size(deltas,2));
    for i = 1:size(deltas,2)
        g_vals(i) = new_tot'*(diag(deltas(i).pr(:,1))*new_tot+deltas(i).pr(:,2));
    end
    [~,m] = max(g_vals); %finds maximum g
    A_deltas = diag(deltas(m).pr(:,1));
    b_deltas = deltas(m).pr(:,2);
    
    
    
    f_cost = x_hat'*(A_0*new_tot+b_0) + new_tot'*(A_deltas * new_tot + b_deltas)/coord.num_ag; 
end