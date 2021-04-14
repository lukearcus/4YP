function[vec_cost] = calc_f_cost(sol_vec,n_ag,T,pr_nom)

vec_cost = zeros(n_ag,1);

x = sol_vec(1:n_ag*T);

x = reshape(x,T,n_ag);
tot_demand = sum(x,2);

% computes common term
for i = 1:n_ag
    vec_cost(i) = x(:,i)'*(pr_nom .* tot_demand);
end

end