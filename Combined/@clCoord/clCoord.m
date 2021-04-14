classdef clCoord < handle
    properties
        tot_demand;
        A_0;
        num_deltas;
        num_ag; % only the minimizing players
        deltamax; % index of current maximizing scenario in the set of samples
        y_max; % linear combination of samples: its final value should be one of the vertices of the simplex
        tot_demand_hat;
        costmax; % stores the max value
        poly;
        maximiser;
    end
    
    methods
        function coord = clCoord(T,pr_nom,n_deltas,n_players,maximise_fun) % constructor
            coord.maximiser = maximise_fun;
            
            coord.tot_demand = zeros(T,1);
            coord.poly = POLYTOPE;
            coord.poly.lb = zeros(n_deltas,1);
            coord.poly.ub = ones(n_deltas,1);
            coord.poly.Aeq = ones(1,n_deltas); %should this be relaxed to inequality constraint?
            coord.poly.beq = 1;
            coord.poly.A = [];
            coord.poly.b = [];
            
            coord.num_deltas = n_deltas;
            coord.num_ag = n_players;

            coord.y_max = zeros(n_deltas,1);
            coord.y_max(1) = 1;
            coord.deltamax = 1; % initialises with the first uncertainty sample

            coord.A_0 = pr_nom; % nominal congestion price
            
            coord.costmax = zeros(n_players,1);
            coord.tot_demand_hat = coord.tot_demand;
        end
        
        function comp_demand(coord,ev)
            
            [coord.tot_demand] = aggr(ev); % computes the aggregate demand according to the updated schedules;
            
        end
        
        function comp_demand_hat(coord,ev)
           [coord.tot_demand_hat] = aggr_hat(ev); 
        end
        
        function maxcentral(coord,coord_fun_inputs) % adversarial player(s) in Facchinei et al 2014
            [new_y_max,flag,f_cost] = coord.maximiser(coord,coord_fun_inputs);
            if any(flag == [1 2])
                
                % updates y_max
                coord.y_max = new_y_max;
                coord.costmax = f_cost;
                
            else
                fprintf('MAX QP Solver error: %d\n',flag)
                beep;
                keyboard;
            end
                        
        end
            
    end

end