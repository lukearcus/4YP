function [costs,ev,coord] = MinMax(setup,samples,ev,coord,stopping,tau,samplesRemoval)
solveropt;
costs = zeros(setup.N_ag+1,stopping.n_iter_MAX,samplesRemoval.k+1);
hist_violation_rate = zeros(samplesRemoval.k,1);
samplesRemoved = 0;
coord.y_max = zeros(size(samples.simple,2),1);
coord.y_max(1) = 1;
while samplesRemoved <= samplesRemoval.k
    t=1;
    coeff = 1;
    differ(t) = stopping.tol_out + 1;
    vec_sched = [];
    for i = 1:setup.N_ag
        vec_sched = [vec_sched; ev(i).schedule];
    end
    vec_cost = zeros(setup.N_ag+1,1); % cost vector (for QP)
    vec_cost_tot = zeros(setup.N_ag+1,1); % actual cost (including independent components that stay out of QP cost)
    
    sol_approx = [vec_sched;coord.y_max]; % solution vector
    differ_cost = stopping.tol_out+1;
    
    cost_out = zeros(setup.N_ag+1,stopping.n_iter_MAX);
    sol_out = zeros(setup.N_ag*setup.T+length(samples.simple),stopping.n_iter_MAX);
    
    coord.comp_demand(ev(1:setup.N_ag));
    while t <= stopping.n_iter_MAX && differ_cost > stopping.tol_out
        
        t = t+1;
        
        fprintf('External loop: %d\n',t);
        
        % updates reg_setpoint (centre of regularization) with the solution of the previous inner loop
        reg_setpoint_col = sol_approx;
        reg_setpoint = reshape(reg_setpoint_col(1:setup.N_ag*setup.T),setup.T,setup.N_ag);
        reg_y = reg_setpoint_col(setup.N_ag*setup.T+1:end);
        
        tt = 1; % initialises inner loop index
        
        differ_inn = 1; % initialises norm difference in solutions of subsequent inner iterations
        
        while tt <= stopping.n_iter_inn_MAX && differ_inn > stopping.tol_inn
            [y_val,y_ind] = max(coord.y_max);
            fprintf('-- Loop %d/%d / d_in = %.4e / d_out = %.5e / D %.5e / max(y) = %.3e / %d\n',t,tt,differ_inn,differ(t-1),differ_cost,y_val,y_ind);
            
            % decentralized algorithm (can be run in parallel)
            for i = 1:setup.N_ag
                % each player only gets aggregate information on other players'
                % strategies, based on previous iteration
                ev_fun_inputs.coord = coord;
                ev_fun_inputs.setpoint = reg_setpoint(:,i);
                ev_fun_inputs.deltas = samples.diag;
                ev_fun_inputs.tau = tau;
                ev_fun_inputs.coeff=coeff;
                ev_fun_inputs.deltaT = setup.deltaT;
                ev_fun_inputs.opt_setup = opt_setup;
                ev(i) = ev(i).better(ev_fun_inputs);
            end
            
            hist_ymax = coord.y_max;
            
            % solution of max
            coord_fun_inputs.ev = ev;
            coord_fun_inputs.reg_y = reg_y;
            coord_fun_inputs.tau=tau;
            coord_fun_inputs.coeff=coeff;
            coord_fun_inputs.opt_setup_max = opt_setup_max;
            coord_fun_inputs.deltas = samples.simple;
            coord.maxcentral(coord_fun_inputs); % updates coord.y_max
            
            % Updates players' decision vector
            hist_sched = vec_sched;
            hist_cost = vec_cost;
            vec_sched = [];
            for i = 1:setup.N_ag
                ev(i).schedule = ev(i).x;
                % composes global vector
                vec_sched = [vec_sched; ev(i).schedule];
                vec_cost(i) = ev(i).cost/setup.N_ag;
            end
            vec_cost(setup.N_ag+1) = coord.costmax/setup.N_ag;
            
            coord.comp_demand(ev(1:setup.N_ag)); % updates total demand (sigma)
            
            sol_approx = [vec_sched;coord.y_max];
            
            differ_inn = norm(vec_cost-hist_cost,inf); % cost difference norm from previous (inner) iteration
            
            tt = tt+1;
        end
        
        differ(t) = norm(sol_approx - reg_setpoint_col,inf); % distance from previous regulation centre
        
        
        cost_prev = vec_cost_tot;
        [vec_cost_tot] = calc_cost(sol_approx,setup.N_ag,setup.T,coord.A_0,samples.col);
        
        cost_out(:,t-1) = vec_cost_tot/setup.N_ag;
        sol_out(:,t-1) = sol_approx;
        
        differ_cost = norm(vec_cost_tot - cost_prev)/setup.N_ag; % cost difference norm from previous (outer) iteration
        
        coeff = 1/(1+t); % coeff is an infinite series converging to zero, and infinite sum
        
        %     diff_y = norm(coord.y_max - y_max_prec);
        y_max_prec = coord.y_max;
    end
    costs(:,:,samplesRemoved+1) = cost_out;
    coord_fun_inputs.ev = ev;
    coord_fun_inputs.reg_y = reg_y;
    coord_fun_inputs.tau=tau;
    coord_fun_inputs.coeff=coeff;
    coord_fun_inputs.opt_setup_max = opt_setup_max;
    coord_fun_inputs.deltas = samples.simple;
    samplesRemoval.comp_card(samplesRemoved+1) = sum(coord.y_max > 1e-5);
    
    hist_violation_rate(samplesRemoved+1) = check_violations(samplesRemoval.tests,samples.simple,coord(1),coord_fun_inputs);
    if samplesRemoval.mode == 1
        [samples.simple,reg_y,removed] = maximum_removal(coord,samples.simple,coord.y_max);
    else
        if samplesRemoved == samplesRemoval.k
            samplesRemoved = samplesRemoved+1;
        end
        [samples.simple,reg_y,removed] = all_active_removal(coord,samples.simple,coord.y_max,samplesRemoval.k-samplesRemoved);
    end
    samplesRemoved = samplesRemoved+removed;
    
    N_samples=size(samples.simple,2);
    samples.col = zeros(N_samples*setup.T,2); % column vectors of uncertain price coefficients (same information of deltas struct but put in column)
    
    temp = [samples.simple(1:N_samples).pr]; % produces matrix with all deltas arranged in column
    
    
    samples.col(:,1) = reshape(temp(:,1:2:end),N_samples*setup.T,1); %  a(t)
    samples.col(:,2) = reshape(temp(:,2:2:end),N_samples*setup.T,1); %  b(t)
    
    % creates diagonal matrices in column (used in bestrespregQP.m), only for a(t), then append b as last column
    samples.diag = zeros(N_samples*setup.T,setup.T);
    
    for jj = 1:N_samples
        ind_i = (jj-1)*setup.T+1;
        ind_o = ind_i+setup.T-1;
        samples.diag(ind_i:ind_o,:) = diag(samples.simple(jj).pr(:,1)); % A_w
    end
    samples.diag = [samples.diag samples.col(:,2)];
    
    coord_fun_inputs.ev = ev;
    coord_fun_inputs.reg_y = reg_y;
    coord_fun_inputs.tau=tau;
    coord_fun_inputs.coeff=coeff;
    coord_fun_inputs.opt_setup_max = opt_setup_max;
    coord_fun_inputs.deltas = samples.simple;
    coord.maxcentral(coord_fun_inputs); % updates coord.y_max
    
    
end
end