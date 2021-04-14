function V = recursive_eps(eps,N,beta,d,k)

for i = 1:1000
    eps = (k + d-1+sqrt(2*eps*N*log(((eps*N)^(d-1))/beta)))/N;
end
V = eps;
end