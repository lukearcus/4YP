function [cost,schedule] = Centralised(setup,deltas,ev,coord,sampleRemoval)
samplesRemoved = 0;
cost = zeros(sampleRemoval.k+1,1);
while samplesRemoved <= sampleRemoval.k
sched = sdpvar(setup.T,setup.N_ag,'full');
sigma = sdpvar(setup.T,1);
deltaVals = sdpvar(size(deltas,2),1);

constraints = [sched(:)>=0,sigma == sum(sched,2)];
for i = 1:setup.N_ag
    constraints = [constraints,sched(:,i)<=ev(i).poly.ub'];
    constraints = [constraints,ev(i).poly.A*sched(:,i)<=ev(i).poly.b]; %maybe relax this?
end
for i = 1:size(deltas,2)
    constraints = [constraints,deltaVals(i) == sigma'*(diag(deltas(i).pr(:,1))*sigma+deltas(i).pr(:,2))];
end

objective = sigma'*(diag(coord.A_0)*sigma); %f objective
objective = objective + max(deltaVals);


options = sdpsettings('verbose',2);

sol = optimize(constraints,objective,options);

schedule = value(sched);
cost(samplesRemoved+1) = value(objective);

%% sample removal
%need to figure out if possible to remove more than one at a time?

gvals = value(deltaVals);

[~,m] = max(gvals);
deltas(m) = [];
samplesRemoved = samplesRemoved+1;
end

end