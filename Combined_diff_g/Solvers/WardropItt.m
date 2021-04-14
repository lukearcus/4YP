function [new_schedule,flag,f_cost,x_hat] = WardropItt(ev,ev_fun_inputs)
    coord = ev_fun_inputs.coord;
    poly = ev.poly;
    z=ev_fun_inputs.z;

    b_0 = zeros(size(coord.tot_demand));
    A_0 = diag(coord.A_0);
    
    T = size(A_0,1);
    A_deltas = ev_fun_inputs.deltas(:,1:end-1);
    b_deltas = ev_fun_inputs.deltas(:,end);
    y_max = coord.y_max;
    num_deltas = size(y_max,1);
    yw_Aw = kron(y_max,ones(T,1)).*A_deltas;
    yw_bw = kron(y_max,ones(T,1)).*b_deltas;
    aux1 = kron(ones(1,num_deltas),eye(T));

    f = (A_0+aux1*yw_Aw)*z+b_0+aux1*yw_bw;
    
    opts = optimoptions('linprog');
    opts = optimoptions(opts,'Display','none');

    [new_schedule,f_cost,flag,~] = linprog(f,poly.A,poly.b,poly.Aeq,poly.beq,poly.lb,poly.ub,opts);
    x_hat = new_schedule;

end