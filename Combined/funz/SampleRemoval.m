%% sets up k and violation checker

epsilon = 0.1;
beta = 1e-5;
d = (N_ag+1)*T;
%d=T;
N=N_samples;
%k = floor(epsilon*N-d+1-sqrt(2*epsilon*N*((d-1)*log(epsilon*N)-log(beta))));
k=20;
%expected_violation_rate = calc_expected_epsilon(k,d,N_samples,beta,0.001);

