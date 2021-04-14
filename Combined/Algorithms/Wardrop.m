function [costs,ev,coord] = Wardrop(setup,samples,ev,coord,stopping,samplesRemoval)
solveropt;
costs = zeros(setup.N_ag,stopping.max_itt,samplesRemoval.k+1);
samplesRemoved = 0;
coord.y_max = zeros(size(samples.simple,2),1);
coord.y_max(1) = 1;
diff = 0;
x_hats = zeros(setup.N_ag,setup.T);
while samplesRemoved <= samplesRemoval.k
    coord.comp_demand(ev(1:setup.N_ag));
    z = coord.tot_demand;
    for h=1:stopping.max_itt
        reg_y = coord.y_max;
        for i = 1:setup.N_ag
            ev_fun_inputs.z = z;
            ev_fun_inputs.coord = coord;
            ev(i) = ev(i).better(ev_fun_inputs);
            costs(i,h,samplesRemoved+1) = ev(i).cost/setup.N_ag;
        end
        for i = 1:setup.N_ag
            ev(i).schedule = ev(i).x;
        end
        coord.comp_demand(ev(1:setup.N_ag));
        z = (1-1/h)*z + (1/h)*coord.tot_demand;
        if h > 2
            diff = norm(costs(:,h-2,samplesRemoved+1)-costs(:,h-1,samplesRemoved+1),inf);
            if( diff < stopping.tol) %maybe check something else for convergence?
                break;
            end
        end
        
%         coord_fun_input.deltas = samples.simple;
%         coord_fun_input.ev = ev;
%         coord_fun_input.reg_y = reg_y;
%         coord_fun_input.tau = 10;
%         coord_fun_input.coeff = 1/(1+h);
%         coord_fun_input.opt_setup_max = opt_setup_max;
%         coord.maxcentral(coord_fun_input); %generates a max_y for use in smaple removal
            g_vals = zeros(size(samples.simple,2),1);
            for m = 1:size(samples.simple,2)
                g_vals(m) = coord.tot_demand'*(samples.simple(m).pr(:,1).*coord.tot_demand + samples.simple(m).pr(:,2));
            end
            coord.y_max = double(g_vals ==  max(g_vals)); 

        approx_sol = [];
        for i =1:setup.N_ag
            x_hats(i,:) = (x_hats(i,:)*(h-1)+ev(i).schedule')/h;
            approx_sol = [approx_sol;x_hats(i,:)'];
        end
        approx_sol = [approx_sol;coord.y_max];
        cost = calc_cost(approx_sol,setup.N_ag,setup.T,coord.A_0,samples.col);
        costs(:,h,samplesRemoved+1) = cost(1:setup.N_ag)/setup.N_ag; %might need tweaking?
        fprintf("current diff: %.4e current max in y: %.5e \n",diff,max(coord.y_max));

    end
    
    if samplesRemoval.mode == 1
        [samples.simple,reg_y,removed] = maximum_removal(coord,samples.simple,coord.y_max);
    else
        if samplesRemoved == k
            samplesRemoved = samplesRemoved+1;
        end
        [samples.simple,reg_y,removed] = all_active_removal(coord,samples.simple,coord.y_max,SamplesRemoval.k-samples);
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
    
%     coord_fun_inputs.ev = ev;
%     coord_fun_inputs.reg_y = reg_y;
%     coord_fun_inputs.tau=10;
%     coord_fun_inputs.coeff=1/(1+h);
%     coord_fun_inputs.opt_setup_max = opt_setup_max;
%     coord_fun_inputs.deltas = samples.simple;
%     coord.maxcentral(coord_fun_inputs); % updates coord.y_max
            coord.y_max = zeros(size(samples.simple,2),1);
        coord.y_max(1) = 1;
end
    
end