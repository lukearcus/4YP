% implements regularized Algorithm 3 from Scutari et al 2014 "Real and
% complex monotone communication games"
% --- solves max problem (player N+1) in Facchinei et al. 2014
function [y_max,exitflag,val_f_cost,output] = maxnashgameQP_v3(coord,coord_fun_inputs)

ev = coord_fun_inputs.ev;
reg_y = coord_fun_inputs.reg_y;
deltas = coord_fun_inputs.deltas;
tau = coord_fun_inputs.tau;
coeff = coord_fun_inputs.coeff;
opt_setup = coord_fun_inputs.opt_setup_max;

% reg_y is the regularization setpoint updated in the outer loop
% tau : proximal regularization coefficient

% deltas is a structure whose elements are Tx2 matrices, representing
% a_d(t) and b_d(t)

% tot_demand,deltas,num_players are known by coordinator
tot_demand = coord.tot_demand; % avg_demand is a vector with values along the horizon T

num_samples = size(reg_y,1);


%% defines unit simplex to describe uncertainty set with samples as vertices
% (see Facchinei et al 2014, (3) and below

%%% constraints: y>=0
low_b = zeros(num_samples,1);
upp_b = Inf*ones(num_samples,1);
%%% 

%%% equality constraints (sum(y) = 1)
Ae = ones(1,num_samples);
be = 1;
%%%

%% defines vector of all possible cost function values due to different samples (vertices)
g = zeros(num_samples,1);

for j = 1:num_samples
    for i = 1:coord.num_ag
        g(j) = g(j) + ev(i).schedule' * (deltas(j).pr(:,1).*tot_demand + deltas(j).pr(:,2));
    end
end


H = sparse(coord.num_ag*(tau + coeff)*eye(num_samples));
g_qp = -g - coord.num_ag*tau*reg_y;

[y_max,val_f_cost,exitflag,output] = quadprog(H,g_qp,[],[],Ae,be,low_b,upp_b,[],opt_setup);

end