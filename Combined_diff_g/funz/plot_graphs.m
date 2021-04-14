function plot_graphs(t,cost_out,N_ag,T,ev,coord,deltas,hist_violation_rate,sample_costs)
    figure('name',"cost");
    norms = ones(t-1,1);
    for i = 1:t-1
       norms(i) = norm(cost_out(:,i)-cost_out(:,t-1)); 
    end
    plot(1:t-1,norms);
    title('cost distance from NE')
    xlabel('number of iterations')
    ylabel('global cost distance from final NE')
    
    schedule = 0;
    figure('name',"schedules & price");
    yyaxis left;
    for i = 1:N_ag
        schedule = schedule+ ev(i).schedule;
    end
    plot(1:T,schedule)
    ylabel('aggregate schedule')
    yyaxis right;
   % [y_val,y_ind] = max(coord.y_max);
    plot(1:T,coord.A_0);
    title('schedules (blue) and energy price (red)')
    xlabel('time slot')
    ylabel('energy price £/kWh')
    
    figure('name',"violation rate and global cost with samples removed")
    
    samples_violation = find(hist_violation_rate)-1;
    hist_violation_rate = nonzeros(hist_violation_rate);
    sample_nums = find(sample_costs)-1;
    sample_costs = nonzeros(sample_costs);
    
    subplot(1,2,1);
    plot(samples_violation,hist_violation_rate);
    title('empirical violation rate')
    ylabel('violation rate')
    xlabel('samples removed')
    subplot(1,2,2);
    plot(sample_nums,sample_costs);
    title('global cost')
    xlabel('samples removed')
    ylabel('global cost')
  %  plot(1:T,deltas(y_ind).pr(:,1));

end