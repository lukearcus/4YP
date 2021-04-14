function[violation_rate] = check_violations(tests,deltas,coord,coord_fun_inputs)
threshold = 1e-5;
violation_rate = 0;
try
if tests == 0
    return
end
catch
    
end
coord_fun_inputs.reg_y = [coord_fun_inputs.reg_y;0];
    prev_y = coord.y_max;
    test = deltas;
for i = 1:size(tests,2)
    test(size(deltas,2)+1)=tests(i);
    
    coord_fun_inputs.deltas = test;
    coord.maxcentral(coord_fun_inputs);
    if coord.y_max(size(test,2)) > threshold
        violation_rate = violation_rate + 1/size(tests,2);
    end
    coord.y_max = prev_y;
end
end