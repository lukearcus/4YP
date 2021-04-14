
clear
%% setup variables
addpath funz Solvers Sample_removal Algorithms


Max_ag = 10;
def_case;
Min_ag = 10;
stopping.tol = 1e-4;
stopping.max_itt = 10000;

stopping2 = stopping;
stopping2.tol = 1e-4;
stoppingW = stopping;
stoppingW.tol = 1e-4;

stoppingMM.n_iter_MAX = 10000;
stoppingMM.n_iter_inn_MAX = 50;
stoppingMM.tol_out = 1e-3;
stoppingMM.tol_inn = 5e-10;

mode = "inc_m";

samplesRemoval.k = 0; %max number of samples to remove
samplesRemoval.mode = 1; %1 for max removal, anything else for all active



sampleSize = 500;

if mode == "Diff_m"
    itt = 1:size(samples.simple,2)/sampleSize;
    val = zeros(size(samples.simple,2)/sampleSize,4,Max_ag,samplesRemoval.k+1,3);
    social_optimals = zeros(size(samples.simple,2)/sampleSize,Max_ag,samplesRemoval.k+1,3);
else
    itt = 2:2:size(samples.simple,2);
    val = zeros(size(samples.simple,2),4,Max_ag,samplesRemoval.k+1,3);
    social_optimals = zeros(size(samples.simple,2),Max_ag,samplesRemoval.k+1,3);
end

timeNR = 0;
timeNR2 = 0;
timeMM = 0;
timeCent = 0;
timeW = 0;

%% run
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
        tau = N_ag;
        setup.N_ag = N_ag;
        for i = 1:N_ag
            
            coordNR(i) = clCoord(i,setup.T,coordNR(1).A_0,size(current_samples.simple,2),N_ag,@maxnashgameQP_v3);
            coordNR2(i) = clCoord(i,setup.T,coordNR(1).A_0,size(current_samples.simple,2),N_ag,@maxnashgameQP_v3);
            coordMM(i) = clCoord(i,setup.T,coordNR(1).A_0,size(current_samples.simple,2),N_ag,@maxnashgameQP_v3);
            coordW(i) = clCoord(i,setup.T,coordNR(1).A_0,size(current_samples.simple,2),N_ag,@maxnashgameQP_v3);
            coordCent(i) = clCoord(i,setup.T,coordNR(1).A_0,size(current_samples.simple,2),N_ag,@maxnashgameQP_v3);
            %
            %             coord(i).num_ag = N_ag;
            %             coord(i).costmax = zeros(N_ag,1);
            %             coord(i).num_deltas=m;
            %             coordNR(i).num_ag = N_ag;
            %             coordNR(i).costmax = zeros(N_ag,1);
            %             coordNR(i).num_deltas=m;
            %             coordNR2(i).num_ag = N_ag;
            %             coordNR2(i).costmax = zeros(N_ag,1);
            %             coordNR2(i).num_deltas=m;
            %             coordMM(i).num_ag = N_ag;
            %             coordMM(i).costmax = zeros(N_ag,1);
            %             coordMM(i).num_deltas=m;
            %             coordW(i).num_ag = N_ag;
            %             coordW(i).costmax = zeros(N_ag,1);
            %             coordW(i).num_deltas=m;
        end
        
        tic;
        if timeNR < 3600
            [costsNR,evNR,coordNR] = NoRegret(setup,current_samples,evNR,coordNR,stopping,0,samplesRemoval);
            timeNR = toc;
        end
        tic;
        if timeNR2 < 3600
            [costsNR2,evNR2,coordNR2] = NoRegret(setup,current_samples,evNR2,coordNR2,stopping2,0,samplesRemoval);
            timeNR2 = toc;
        end
        tic;
        if timeMM < 3600
            [costsMM,evMM,coordMM] = MinMax(setup,current_samples,evMM,coordMM,stoppingMM,tau,samplesRemoval);
            timeMM = toc;
        end
        tic;
        if timeCent < 3600
            [costsCent,scheduleCent] = Centralised(setup,current_samples,evNR,coordNR,samplesRemoval);
            timeCent = toc;
        end
        tic;
        if timeW < 3600
            [costsWard,evW,coordW] = Wardrop(setup,current_samples,evW,coordW,stoppingW,samplesRemoval);
            timeW = toc;
        end
        for j = 1:samplesRemoval.k+1
            for k = 1:3
                social_optimum = costsCent(j,k)/N_ag;
                social_optimals(m,N_ag,j,k) = social_optimum;
                if timeNR < 3600
                    val(m,1,N_ag,j,k) = (sum(costsNR(:,nnz(costsNR(1,:,j,k)),j)))/social_optimum;
                else
                    val(m,1,N_ag,j,k) = inf;
                end
                if timeNR2 < 3600
                    val(m,2,N_ag,j,k) = (sum(costsNR2(:,nnz(costsNR2(1,:,j,k)),j)))/social_optimum;
                else
                    val(m,2,N_ag,j,k) = inf;
                end
                if timeMM <3600
                    val(m,3,N_ag,j,k) = (sum(costsMM(1:N_ag,nnz(costsMM(1,:,j,k)),j)))/social_optimum;
                else
                    val(m,3,N_ag,j,k) = inf;
                end
                if timeW < 3600
                    val(m,4,N_ag,j,k) = (sum(costsWard(1:N_ag,nnz(costsWard(1,:,j,k))-1,j)))/social_optimum;
                else
                    val(m,4,N_ag,j,k) = inf;
                end
            end
        end
        
        %         if N_ag+1 < Max_ag
        %             stopping.tol = (stopping.tol*(N_ag+1))/(N_ag);
        %             stopping2.tol = (stopping2.tol*(N_ag+1))/(N_ag);
        %             stoppingMM.tol_out = (stoppingMM.tol_out*(N_ag+1))/(N_ag);
        %             stoppingMM.tol_inn = (stoppingMM.tol_inn*(N_ag+1))/(N_ag);
        %             stoppingW.tol = (stoppingW.tol*(N_ag+1))/(N_ag);
        %         end
    end
end