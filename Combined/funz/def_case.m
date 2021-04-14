% initialises random number generator
rng('default');
rng(2);

% number of players % 
N_ag = Max_ag;

T = 6; % num of time slots composing the optimization interval
deltaT = 4; % duration of a single time slot [h]

rate_limit = randi([6,15],N_ag,1); % [kW] max charge rate by each user

energy_amount = 0.05*T*deltaT*rate_limit.*rand(N_ag,1); % [kWh] 

efficiency = 0.8+0.2*rand(N_ag,1); %minimum efficiency of 0.8

constraints = ones(N_ag,T);
constraints = constraints.*rate_limit;

starts = randi([1,T],N_ag,1);
lens = randi([0,round(T/5)],N_ag,1);



teta_in = [];
teta_len = [];

% assegna parametri a EV
for i = 1:N_ag
    constraint = constraints(i,:);
    constraint(starts(i):min(starts(i)+lens,T)) = 0;
    evNR(i) = clPlayer(i,T,constraint,energy_amount(i),deltaT,efficiency(i),@SubGradientDescent);
    evNRfixed(i) = clPlayer(i,T,constraint,energy_amount(i),deltaT,efficiency(i),@SubGradientDescent);
    evNRGrad(i) = clPlayer(i,T,constraint,energy_amount(i),deltaT,efficiency(i),@GradientDescent);
    evNRGradfixed(i) = clPlayer(i,T,constraint,energy_amount(i),deltaT,efficiency(i),@GradientDescent);
    evMM(i) = clPlayer(i,T,constraint,energy_amount(i),deltaT,efficiency(i),@bestrespregQP_v3);
    evW(i) = clPlayer(i,T,constraint,energy_amount(i),deltaT,efficiency(i),@WardropItt);
end

%% energy price

energy_price = [ % gio 19 gen 2017 (from bmreports, 48 half-hour slots) thu 19 jan 2017
43.81
44.2
43.68
43.43
42.3
42.71
42.07
41.92
41.72
41.28
42
42.34
47.69
42.73
51.38
59.59
62.06
57.55
55.17
56.33
55.38
53.1
46.38
45.55
44.98
44.64
43.93
43.43
43.17
42.45
49.78
51.95
53.3
80.72
69.88
68.43
71.23
67.64
57.58
52.13
50.31
49.13
49.65
48.17
47.58
46.47
43.1
42.67]; % [£/MWh]

energy_price = energy_price/1e3; % converts to [£/kWh]

%% uncertainty samples
% uncertainty is on affine pricing function coefficients: a_d(t)*sigma(x) + b_d
% samples are in the form {a_j(t),b_j(t)}
% arranges samples in data structures used later in the code

N_samples = 500; % n samples

% nominal b(t) from daily price profile; nominal a(t) = constant (defined
% in the coordinator object constructor)
A_0 = energy_price(1:48/T:end); % b_0 is assumed = 0


for j = 1:N_samples
    deltas(j).pr = [2e-2*exp(randn(T,1))*20/N_ag    2e-2*(rand(T,1)-0.5)];
end

deltas_col = zeros(N_samples*T,2); % column vectors of uncertain price coefficients (same information of deltas struct but put in column)

temp = [deltas(1:N_samples).pr]; % produces matrix with all deltas arranged in column

deltas_col(:,1) = reshape(temp(:,1:2:end),N_samples*T,1); %  a(t)
deltas_col(:,2) = reshape(temp(:,2:2:end),N_samples*T,1); %  b(t)

% creates diagonal matrices in column (used in bestrespregQP.m), only for a(t), then append b as last column
deltas_diag = zeros(N_samples*T,T); 

for jj = 1:N_samples
    ind_i = (jj-1)*T+1;
    ind_o = ind_i+T-1;
    deltas_diag(ind_i:ind_o,:) = diag(deltas(jj).pr(:,1)); % A_w
end
deltas_diag = [deltas_diag deltas_col(:,2)];

%% initializes the coordinator
coordNR = clCoord(T,A_0,N_samples,N_ag,@maxnashgameQP_v3);
coordNRfixed = clCoord(T,A_0,N_samples,N_ag,@maxnashgameQP_v3);
coordNRGrad= clCoord(T,A_0,N_samples,N_ag,@maxnashgameQP_v3);
coordNRGradfixed= clCoord(T,A_0,N_samples,N_ag,@maxnashgameQP_v3);
coordMM = clCoord(T,A_0,N_samples,N_ag,@maxnashgameQP_v3);
coordW = clCoord(T,A_0,N_samples,N_ag,@maxnashgameQP_v3);

%% sample removal stuff

epsilon = 0.1;
beta = 1e-5;
%d = (N_ag+1)*T;
d=T;
N=size(deltas,2);
%k = floor(epsilon*N-d+1-sqrt(2*epsilon*N*((d-1)*log(epsilon*N)-log(beta))));

samplesRemoval.k = 0; %max number of samples to remove
N_tests = 000; %number of tests for checking violation rate
samplesRemoval.mode = 1; %1 for max removal, anything else for all active
%samplesRemoval.expectedEps = calc_expected_epsilon(samplesRemoval.k,d,N,beta,0.01);
samplesRemoval.comp_card = zeros(1,samplesRemoval.k+1);
for j = 1:N_tests
    tests(j).pr = [2e-2*exp(randn(T,1))*20/N_ag    2e-2*(rand(T,1)-0.5)];
end
if N_tests > 0
    samplesRemoval.tests = tests;
else
    samplesRemoval.tests=0;
end
%% cleans up
Setup.A_0 = A_0;
setup.deltaT = deltaT;
setup.N_ag = N_ag;
setup.T = T;

samples.simple = deltas;
samples.col = deltas_col;
samples.diag = deltas_diag;
clearvars -except coordNR evNR coordNRfixed evNRfixed coordNRGrad evNRGrad coordNRGradfixed evNRGradfixed coordMM evMM coordW evW samples setup Max_ag samplesRemoval