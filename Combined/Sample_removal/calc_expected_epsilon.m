function epsilon = calc_expected_epsilon(k,d,N,beta,accuracy)
%using sampling and dsicarding approach, epsilon produced is directly
%violation rate
beta_calc = 1;
epsilon = 0;
warning('off','MATLAB:nchoosek:LargeCoefficient')
while beta_calc > beta
    epsilon = epsilon+accuracy;
    beta_calc = 0;
        for i = 0:k+d-1
           beta_calc =  beta_calc + nchoosek(N,i)*(epsilon^i)*(1-epsilon)^(N-i);
        end
        beta_calc = beta_calc*nchoosek(k+d-1,k);
end
warning('on','MATLAB:nchoosek:LargeCoefficient')

end