function [new_schedule,flag,f_cost,x_hat,eta_sum] = SubGradientDescent(ev,ev_fun_inputs)
    coord = ev_fun_inputs.coord;
    deltas = ev_fun_inputs.deltas;
    poly = ev.poly;
    stepsize = ev_fun_inputs.coeff;
    t = ev_fun_inputs.t;
    tot_demand = coord.tot_demand;
    x_meno_i = tot_demand - ev.schedule;
    b_0 = zeros(size(tot_demand));
    curr_sched = ev.schedule;
    A_0 = diag(coord.A_0);
    num_ag = coord.num_ag;
    f_gradient = 2*A_0*ev.schedule + A_0*x_meno_i+b_0;
      

    g_vals = zeros(size(deltas,2),1);
    for i = 1:size(deltas,2)
        g_vals(i) = tot_demand'*(diag(deltas(i).pr(:,1))*tot_demand+deltas(i).pr(:,2));
    end
    [~,m] = max(g_vals); %finds maximum g
    A_deltas = diag(deltas(m).pr(:,1)); %uses A corresponding to max g (effectively uses a subgradient)
    b_deltas = deltas(m).pr(:,2);
    g_gradient = 2*A_deltas*tot_demand+b_deltas;
    
    gradient = f_gradient+g_gradient/num_ag;
    sched = curr_sched - stepsize*gradient;
    
    new_schedule = poly.project_onto(sched);
    flag = 1;
    
    %cost should be defined using x_hat
    x_hat = (ev.x_hat*(t-1) + new_schedule)/t;
    %x_hat = ((ev.x_hat*ev.eta_sum) + new_schedule*stepsize)/(ev.eta_sum + stepsize);
    %x_hat = poly.project_onto(x_hat); %to ensure x_hat remains in feasible set
    eta_sum = ev.eta_sum + stepsize;
    
    new_tot = coord.tot_demand_hat - ev.x_hat + x_hat;
    
    for i = 1:size(deltas,2)
        g_vals(i) = new_tot'*(diag(deltas(i).pr(:,1))*new_tot+deltas(i).pr(:,2));
    end
    
    f_cost = x_hat'*(A_0*new_tot+b_0) + max(g_vals)/num_ag; 
end