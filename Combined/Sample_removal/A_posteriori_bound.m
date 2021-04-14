function epsilon = A_posteriori_bound(M,k,beta)
warning('off','MATLAB:nchoosek:LargeCoefficient')
epsilon = 1-(beta/(M*nchoosek(M,k)))^(1/(M-k));
warning('on','MATLAB:nchoosek:LargeCoefficient')
end