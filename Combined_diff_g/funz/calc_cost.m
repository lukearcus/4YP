function [vec_cost] = calc_cost(sol_vec,n_ag,T,pr_nom,deltas_col)

vec_cost = zeros(2*n_ag,1);

x = sol_vec(1:n_ag*T);
y = sol_vec(n_ag*T+1:end);

num_deltas = size(y,1)/n_ag;

a_deltas = deltas_col(:,1);
b_deltas = deltas_col(:,2);

a_deltas = reshape(a_deltas,T,num_deltas);
b_deltas = reshape(b_deltas,T,num_deltas);

x = reshape(x,T,n_ag);
y = reshape(y,num_deltas,n_ag);
tot_demand = sum(x,2);

% computes common term
g = zeros(num_deltas,1);
g_cost = zeros(n_ag,1);
for i = 1:n_ag
    for j = 1:num_deltas
        g(j) = x(:,i)'*(a_deltas(:,j) .* tot_demand + b_deltas(:,j));
        
    end
    g_cost(i) = max(g);
    vec_cost(n_ag + i) = y(:,i)'*g;
end


for i = 1:n_ag
    vec_cost(i) = x(:,i)'*(pr_nom .* tot_demand) + g_cost(i);
end