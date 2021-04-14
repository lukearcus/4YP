clear
addpath funz Solvers Sample_removal Algorithms
Max_ag = 10;
def_case;
stopping.tol = 1e-6;
stopping.max_itt = 10000;

stoppingfixed = stopping;
stoppingfixed.tol = 1e-6;
stoppingW = stopping;
stoppingW.tol = 5e-5;

% val = zeros(6,Max_ag,samplesRemoval.k+1);

tau = 10;
stoppingMM.n_iter_MAX = 40000;
stoppingMM.n_iter_inn_MAX = 500;
stoppingMM.tol_out = 1e-6;
stoppingMM.tol_inn = 5e-10;

Min_ag=10;

samplesRemoval.mode = 1;

mode = "inc_m";
sampleSize = 500;

if mode == "Diff_m"
   itt = 1:size(samples.simple,2)/sampleSize;
   val = zeros(size(samples.simple,2)/sampleSize,4,Max_ag,samplesRemoval.k+1,3);
   social_optimals = zeros(size(samples.simple,2)/sampleSize,Max_ag,samplesRemoval.k+1,3);
else
    itt = 2:2:size(samples.simple,2);
    val = zeros(size(samples.simple,2),6,Max_ag,samplesRemoval.k+1);
   social_optimals = zeros(size(samples.simple,2),Max_ag,samplesRemoval.k+1);
end

for m = itt
    if mode == "Diff_m"
        %for different m
        current_samples.simple = samples.simple((m-1)*sampleSize+1:m*sampleSize);
        current_samples.col = samples.col((m-1)*sampleSize*setup.T+1:setup.T*m*sampleSize,:);
        current_samples.diag = samples.diag((m-1)*sampleSize*setup.T+1:setup.T*m*sampleSize,:);
    else
        %increasing m
        current_samples.simple = samples.simple(1:m);
        current_samples.col = samples.col(1:m*setup.T,:);
        current_samples.diag = samples.diag(1:m*setup.T,:);
    end
    for N_ag = Min_ag:Max_ag
        setup.N_ag = N_ag;
        coord.num_ag = N_ag;
        coord.costmax = zeros(N_ag,1);
        coordNR.num_ag = N_ag;
        coordNR.costmax = zeros(N_ag,1);
        coordNRfixed.num_ag = N_ag;
        coordNRfixed.costmax = zeros(N_ag,1);
        coordNRGrad.num_ag = N_ag;
        coordNRGrad.costmax = zeros(N_ag,1);
        coordNRGradfixed.num_ag = N_ag;
        coordNRGradfixed.costmax = zeros(N_ag,1);
        coordMM.num_ag = N_ag;
        coordMM.costmax = zeros(N_ag,1);
        coordW.num_ag = N_ag;
        coordW.costmax = zeros(N_ag,1);
        
        tic;
        [costsNR,evNR,coordNR] = NoRegret(setup,current_samples,evNR,coordNR,stopping,0,samplesRemoval);
        NRtime = toc;
        tic;
        %     [costsNRfixed,evNRfixed,coordNRfixed] = NoRegret(setup,current_samples,evNRfixed,coordNRfixed,stoppingfixed,1,samplesRemoval);
        NRfixedtime = toc;
        tic;
        [costsNRGrad,evNRGrad,coordNRGrad] = NoRegret(setup,current_samples,evNRGrad,coordNRGrad,stopping,0,samplesRemoval);
        NRGradtime = toc;
        tic;
        %     [costsNRGradfixed,evNRGradfixed,coordNRGradfixed] = NoRegret(setup,current_samples,evNRGradfixed,coordNRGradfixed,stoppingfixed,1,samplesRemoval);
        NRGradfixedtime = toc;
        tic;
        [costsMM,evMM,coordMM] = MinMax(setup,current_samples,evMM,coordMM,stoppingMM,tau,samplesRemoval);
        MMtime = toc;
        
        [costsCent,scheduleCent] = Centralised(setup,current_samples.simple,evNR,coordNR,samplesRemoval);
        tic;
        [costsWard,evW,coordW] = Wardrop(setup,current_samples,evW,coordW,stoppingW,samplesRemoval);
        wardTime = toc;
        for j = 1:samplesRemoval.k+1
            social_optimals(m,N_ag,j) = costsCent(j)/N_ag;
            val(m,1,N_ag,j) = (sum(costsNR(:,nnz(costsNR(1,:,j)),j)))/social_optimals(m,N_ag,j);
            %         val(2,N_ag,j) = (sum(costsNRfixed(:,nnz(costsNRfixed(1,:,j)),j)))/social_optimals(m,N_ag,j);
            val(m,3,N_ag,j) = (sum(costsNRGrad(:,nnz(costsNRGrad(1,:,j)),j)))/social_optimals(m,N_ag,j);
            %val(4,N_ag,j) = (sum(costsNRGradfixed(:,nnz(costsNRGradfixed(1,:,j)),j)))/social_optimals(m,N_ag,j);
            val(m,5,N_ag,j) = (sum(costsMM(1:N_ag,nnz(costsMM(1,:,j)),j)))/social_optimals(m,N_ag,j);
            val(m,6,N_ag,j) = (sum(costsWard(1:N_ag,nnz(costsWard(1,:,j))-1,j)))/social_optimals(m,N_ag,j);
        end
    end
end