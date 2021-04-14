function[violation_rate] = check_violations(tests,deltas,coord,coord_inputs,N_ag)
violation_rate = 0;
try
if tests == 0
    return
end
catch
end
threshold = 1e-4;
prev_y = zeros(size(coord(1).y_max,1),N_ag);
for i = 1:N_ag
    prev_y(:,i) = coord(i).y_max;
end
test = deltas;
coord_fun_inputs = coord_inputs;
reg_y = coord_inputs.reg_y;
for i = 1:size(tests,2)
    test(size(deltas,2)+1)=tests(i);
    
    coord_fun_inputs.deltas = test;
   
    for j = 1:N_ag
        coord_fun_inputs.reg_y = [reg_y(:,j);0];
        coord(j).maxcentral(coord_fun_inputs);
        if coord(j).y_max(size(test,2)) > threshold
            violation_rate = violation_rate + 1/size(tests,2);
            coord(j).y_max = prev_y(:,j);
            break;
        end
        coord(j).y_max = prev_y(:,j);
    end
end
end