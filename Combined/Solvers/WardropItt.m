function [new_schedule,flag,f_cost,x_hat] = WardropItt(ev,ev_fun_inputs)
    coord = ev_fun_inputs.coord;
    poly = ev.poly;
    z=ev_fun_inputs.z;

    b_0 = zeros(size(coord.tot_demand));
    A_0 = diag(coord.A_0);

    f = A_0*z+b_0;
    
    opts = optimoptions('linprog');
    opts = optimoptions(opts,'Display','none');

    [new_schedule,f_cost,flag,~] = linprog(f,poly.A,poly.b,poly.Aeq,poly.beq,poly.lb,poly.ub,opts);
    x_hat = new_schedule;

end