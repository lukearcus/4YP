function [removed,deltas,coord,violation] = sampleRemoval(deltas,mode,coord,N_tests,ev)
violation = check_violations(N_tests,deltas,coord,ev,reg_y,tau,coeff,opt_setup_max);

end