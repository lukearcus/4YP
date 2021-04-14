classdef clPlayer < handle
    
    properties
        id;
        schedule ; % schedule over the period
        x; % temp schedule
        cost; % value of the cost function at the last iteration
        poly;
        grad;
        x_hat;
        eta_sum;
        minimiser;

    end
    
	methods
        function ev = clPlayer(id,T,constraints,en_amount,deltaT,efficiency,minimiser_fun) % constructor
            ev.minimiser = minimiser_fun;
            ev.id = id; % id of the player
            ev.poly = POLYTOPE;
            ev.poly.lb = zeros(T,1);
            ev.poly.ub = constraints;
            ev.poly.Aeq = []; %should this be relaxed to inequality constraint?
            ev.poly.beq = [];
            ev.poly.A = -ones(1,T)*efficiency*deltaT;
            ev.poly.b = -en_amount;
          
            ev.schedule = zeros(T,1); % this is updated only at the end of the inner loop (to implement Jacobi algorithm)

            % initialises with a feasible solution
            ev.schedule = ones(T,1) * en_amount/(efficiency*T*deltaT); % deltaT [h] is time slot length 
            ev.schedule = ev.poly.project_onto(ev.schedule);
            ev.cost = 0; %latest value of cost function
            ev.x_hat = ev.schedule;
            
            ev.x = ev.schedule; % new schedule is temporarily stored here (until the end of the inner loop)
            
        end
        
        
        function ev = better(ev,ev_fun_inputs)
            %coord,deltas,coeff
            % This implements Algorithm 3 in Scutari et al 2014 "Real and
            % complex monotone communication games"
            [new_schedule,flag,f_cost,ev.x_hat] = ev.minimiser(ev,ev_fun_inputs);
            
            %coord,deltas,ev.poly,coeff
            %[new_schedule,flag,f_cost,ev.grad,ev.x_hat,ev.eta_sum] = GradientDescent(ev,coord.tot_demand,coord.A_0,coord.num_ag,...
                %deltas,ev.poly,t,coord.tot_demand_hat);
            
            if any(flag == [1 2])  % for meaning of flags check Matlab help (in this case for quadprog)
                
                % updates the schedule
                    ev.x = new_schedule;
                    ev.cost = f_cost;
            else
                fprintf('Min QP Solver error: %d\n',flag)
                beep;
                keyboard
            end
                        
        end
    end
    
end
        
            