function [costs,ev,coord] = NoRegret(setup,deltas,ev,coord,stopping,mode,samplesRemoval)
costs = zeros(setup.N_ag,stopping.max_itt,samplesRemoval.k+1,3);
hist_violation_rate = zeros(1,samplesRemoval.k);
samplesRemoved = 0;
for i = 1:setup.N_ag
    coord(i).y_max = zeros(size(deltas.simple,2),1);
    coord(i).y_max(1) = 1;
end
diff = 10;
while samplesRemoved <= samplesRemoval.k
    %grads = zeros(setup.N_ag,stopping.max_itt);
    for i = 1:setup.N_ag
        coord(i).comp_demand(ev);
    end
    t=0;
    if mode ~=0
        coeff = calc_step_size(ev(1:setup.N_ag),deltas.simple,coord)*0.9; %calc_step_size gives max stepsize, reduce to ensure within bound
        if coeff == 0
           coeff = 1; 
        end
        %coeff is stepsize for descent
    end
    tau = 10;
    solveropt;
    for i = 1:setup.N_ag
        coord(i).comp_demand(ev(1:setup.N_ag));
        coord(i).comp_demand_hat(ev(1:setup.N_ag));
    end
    for itt = 1:stopping.max_itt
        t = t+1;
        if mode == 0
            coeff = 1/sqrt(t);
        end
        
        for i = 1:setup.N_ag
            reg_y(i,:) = coord(i).y_max;
        end
        for i = 1:setup.N_ag
            ev_fun_inputs.coord = coord(i);
            ev_fun_inputs.deltas = deltas.simple;
            ev_fun_inputs.coeff = coeff;
            ev_fun_inputs.t = t;
            ev(i) = ev(i).better(ev_fun_inputs);
            costs(i,t,samplesRemoved+1,1) = ev(i).cost/setup.N_ag;
        end
        %     ev(1).x = [0 100 0 0 0 0]'; %makes agent 1 lie about its preferred
        %     charging strategy
        %     ev(1).x_hat = ev(1).x;
        sol_approx = [];
        for i = 1:setup.N_ag
            ev(i).schedule = ev(i).x;
            sol_approx = [sol_approx;ev(i).x];
        end
        for i = 1:setup.N_ag
            coord(i).comp_demand(ev(1:setup.N_ag));
            coord(i).comp_demand_hat(ev(1:setup.N_ag));
        end
        if t > 1
            diff = norm(costs(:,t-1,samplesRemoved+1)-costs(:,t,samplesRemoved+1),inf);
            if(diff < stopping.tol)
                break;
            end
        end
        

        
        if isequal(ev(1).minimiser,@GradientDescent)
            for i = 1:setup.N_ag
                coord_fun_input.deltas = deltas.simple;
                coord_fun_input.ev = ev;
                coord_fun_input.reg_y = reg_y(i,:)';
                coord_fun_input.tau = tau;
                coord_fun_input.coeff = 1/(1+t);
                coord_fun_input.opt_setup_max = opt_setup_max;
                coord(i).maxcentral(coord_fun_input);
            end
        else
            for i = 1:setup.N_ag
            g_vals = zeros(size(deltas.simple,2),1);
            for m = 1:size(deltas.simple,2)
                g_vals(m) = ev(i).schedule'*(deltas.simple(m).pr(:,1).*coord(i).tot_demand + deltas.simple(m).pr(:,2));
            end
            coord(i).y_max = double(g_vals ==  max(g_vals));
            end
        end
        
        fprintf("current diff: %.4e\n",diff);
        
        
        costs(:,t,samplesRemoved+1,2)=calc_alt_cost(sol_approx,setup.N_ag,setup.T,coord(1).A_0);
        costs(:,t,samplesRemoved+1,3)=calc_f_cost(sol_approx,setup.N_ag,setup.T,coord(1).A_0);

    end
    
    %% code below is for dishonest agent
    %
    % ev(1).schedule = ev(1).poly.project_onto(ones(6,1) * ev(1).poly.b/ev(1).poly.b); %replace 6 with number of timeslots
    % ev(1).x_hat = ev(1).schedule;
    % coord.comp_demand(ev(1:setup.N_ag));
    % coord.tot_demand_hat = coord.tot_demand;
    %
    % for k = 1:100
    %         t = t+1;
    %         coeff = 1/sqrt(t);
    %         ev_fun_inputs.coord = coord;
    %         ev_fun_inputs.deltas = deltas.simple;
    %         ev_fun_inputs.coeff = coeff;
    %         ev_fun_inputs.t = t;
    %         ev(1) = ev(1).better(ev_fun_inputs);
    %         ev(1).schedule = ev(1).x;
    %
    %         coord.comp_demand(ev(1:setup.N_ag));
    %         coord.comp_demand_hat(ev(1:setup.N_ag));
    % end
    %
    % coord.comp_demand(ev(1:setup.N_ag));
    % g = zeros(1,size(deltas.simple,2));
    % for j = 1:size(deltas.simple,2)
    %     g(j) = coord.tot_demand'*(deltas.simple(j).pr(:,1).*coord.tot_demand+deltas.simple(j).pr(:,2));
    % end
    % for i = 1:coord.num_ag
    %     costs(i,200) = ev(i).schedule'*(coord.A_0.*coord.tot_demand) + max(g)/coord.num_ag;
    % end
    
    %% sample removal
    
    %hist_violation_rate(samplesRemoved+1) = check_violations(samplesRemoval.N_tests,samples.simple,coord,ev,reg_y,tau,1,opt_setup_max); %this might need to be tweaked?
    %change violation rate check
    y_vals = zeros(size(deltas.simple,2),1);
    for i = 1:setup.N_ag
       y_vals = y_vals + coord(i).y_max; 
    end
    y_vals = y_vals./setup.N_ag;
    if samplesRemoval.mode == 1
        [deltas.simple,reg_y,removed] = maximum_removal(coord,deltas.simple,y_vals,reg_y);
    else
        if samplesRemoved == k
            samplesRemoved = samplesRemoved+1;
        end
        [deltas.simple,reg_y,removed] = all_active_removal(coord(1),deltas.simple,coord(1).y_max,SamplesRemoval.k-samples);
    end
    samplesRemoved = samplesRemoved+removed;
    
    N_samples=size(deltas.simple,2);
    deltas.col = zeros(N_samples*setup.T,2); % column vectors of uncertain price coefficients (same information of deltas struct but put in column)
    
    temp = [deltas.simple(1:N_samples).pr]; % produces matrix with all deltas arranged in column
    
    
    deltas.col(:,1) = reshape(temp(:,1:2:end),N_samples*setup.T,1); %  a(t)
    deltas.col(:,2) = reshape(temp(:,2:2:end),N_samples*setup.T,1); %  b(t)
    
    % creates diagonal matrices in column (used in bestrespregQP.m), only for a(t), then append b as last column
    deltas.diag = zeros(N_samples*setup.T,setup.T);
    
    for jj = 1:N_samples
        ind_i = (jj-1)*setup.T+1;
        ind_o = ind_i+setup.T-1;
        deltas.diag(ind_i:ind_o,:) = diag(deltas.simple(jj).pr(:,1)); % A_w
    end
    deltas.diag = [deltas.diag deltas.col(:,2)];
    
    for i = 1:setup.N_ag
    if isequal(ev(1).minimiser,@GradientDescent)
        coord_fun_inputs.ev = ev;
        coord_fun_inputs.reg_y = reg_y(i,:)';
        coord_fun_inputs.tau=tau;
        coord_fun_inputs.coeff=1/(1+t);
        coord_fun_inputs.opt_setup_max = opt_setup_max;
        coord_fun_inputs.deltas = deltas.simple;
        coord(i).maxcentral(coord_fun_inputs); % updates coord.y_max
    else
        coord(i).y_max = zeros(size(deltas.simple,2),1);
        coord(i).y_max(1) = 1;
    end
    end
end


end