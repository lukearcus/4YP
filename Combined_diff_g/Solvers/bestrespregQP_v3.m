% implements best response for individual player
function [new_schedule,exitflag,val_f_cost,x_hat,eta_sum] = bestrespregQP_v3(ev,ev_fun_inputs)
% pr_nom is a column vector
coord = ev_fun_inputs.coord;
A_0 = coord.A_0;
num_ag = coord.num_ag;
reg_setpoint = ev_fun_inputs.setpoint;
y_max = coord.y_max;
tau = ev_fun_inputs.tau;
coeff = ev_fun_inputs.coeff;
deltaT = ev_fun_inputs.deltaT;
opt_setup = ev_fun_inputs.opt_setup;
A_deltas = ev_fun_inputs.deltas(:,1:end-1);
b_deltas = ev_fun_inputs.deltas(:,end);

tot_demand = coord.tot_demand;
% y_max is the variable maximised by the adversarial player N+1 (Facchinei et al 2014)

% tau: proximal regularization step

T = size(tot_demand,1);
num_deltas = size(y_max,1); % number of samples

% computes total demand of other players ONLY (removes component of
% current player from the total demand)
tot_meno_i = tot_demand - ev.schedule; 

yw_Aw = kron(y_max,ones(T,1)).*A_deltas;
yw_bw = kron(y_max,ones(T,1)).*b_deltas;

aux1 = kron(ones(1,num_deltas),eye(T)); % sums matrices over w

H = diag(A_0) + aux1*yw_Aw;
h1 = A_0 .* tot_meno_i;
h2 = aux1*(yw_Aw * tot_meno_i + yw_bw);

% constraints
low_b = ev.poly.lb;
upp_b = ev.poly.ub;

% equality constraint relaxed to (active) inequality one 
A = ev.poly.A; % constraint on total energy demand
b = ev.poly.b;

[new_schedule,val_f_cost,exitflag] = quadprog(2*H + tau*eye(T) + coeff*eye(T),h1+h2-tau*reg_setpoint,A,b,[],[],low_b,upp_b,[],opt_setup);

x_hat = new_schedule;
eta_sum = 1; %change this?

end