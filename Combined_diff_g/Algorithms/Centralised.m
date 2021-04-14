function [cost,schedule] = Centralised(setup,samples,ev,coord,sampleRemoval)
samplesRemoved = 0;
cost = zeros(sampleRemoval.k+1,3);
deltas = samples.simple;
while samplesRemoved <= sampleRemoval.k
    sched = sdpvar(setup.T,setup.N_ag,'full');
    sigma = sdpvar(setup.T,1);
    deltaVals = sdpvar(size(deltas,2),setup.N_ag,'full');
    
    constraints = [sched(:)>=0,sigma == sum(sched,2)];
    for i = 1:setup.N_ag
        constraints = [constraints,sched(:,i)<=ev(i).poly.ub'];
        constraints = [constraints,ev(i).poly.A*sched(:,i)<=ev(i).poly.b];
        for j = 1:size(deltas,2)
            constraints = [constraints,deltaVals(j,i) == sched(:,i)'*(deltas(j).pr(:,1).*sigma+deltas(j).pr(:,2))];
        end
    end
    
    objective = sigma'*(diag(coord(1).A_0)*sigma); %f objective
    for i = 1:setup.N_ag
        objective = objective + max(deltaVals(:,i));
    end
    
    
    options = sdpsettings('verbose',2);
    
    sol = optimize(constraints,objective,options);
    
    schedule = value(sched);
    sol_approx = reshape(schedule,[setup.T*setup.N_ag,1]);
    cost(samplesRemoved+1,1) = value(objective);
    cost(samplesRemoved+1,2)=sum(calc_alt_cost(sol_approx,setup.N_ag,setup.T,coord(1).A_0));
    cost(samplesRemoved+1,3)=sum(calc_f_cost(sol_approx,setup.N_ag,setup.T,coord(1).A_0));

    %% sample removal
    %need to figure out if possible to remove more than one at a time?
    
    gvals = value(deltaVals);
    y_vals = zeros(size(deltas,2),1);
    for i =1:setup.N_ag
        y_vals = y_vals+gvals(:,i)';
    end
    y_vals = y_vals./setup.N_ag;
    [~,m] = max(y_vals);
    deltas(m) = [];
    samplesRemoved = samplesRemoved+1;
end

end