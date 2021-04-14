function [vec_cost] = calc_alt_cost(sol_vec,n_ag,T,pr_nom)

x = sol_vec(1:n_ag*T);
% y = sol_vec(n_ag*T+1:end);

% num_deltas = size(deltas_col,1)/T;
vec_cost = zeros(n_ag,1);

% a_deltas = deltas_col(:,1);
% b_deltas = deltas_col(:,2);
% 
% a_deltas = reshape(a_deltas,T,num_deltas);
% b_deltas = reshape(b_deltas,T,num_deltas);

x = reshape(x,T,n_ag);
% if size(y,1) >0
%     y = reshape(y,num_deltas,n_ag);
% end
tot_demand = sum(x,2);

% computes common term
% g = zeros(num_deltas,1);
% g_cost = zeros(n_ag,1);
% for i = 1:n_ag
%     for j = 1:num_deltas
%         g(j) = x(:,i)'*(a_deltas(:,j) .* tot_demand + b_deltas(:,j));
%         
%     end
%     g_cost(i) = max(g);
%     if size(y,1) > 0
%         vec_cost(n_ag + i) = y(:,i)'*g;
%     end
% end

rho = [ 0 0 0 0 2.4 0]';

over_use = tot_demand > rho;

k = 3;

for i = 1:n_ag
    vec_cost(i) = x(:,i)'*(pr_nom .* tot_demand) + sum(over_use*k)/n_ag;
end