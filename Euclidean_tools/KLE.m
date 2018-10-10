function [U,S,V,Prnd,Pmu,f,z] = KLE(P,nrnd)

Pmu = mean(P,1)';
N = size(P,1); m = size(P,2);

%% build orthonormal basis U for KL expansion
[U,S,V] = svd(1/sqrt(N)*P', 'econ');

%% empirical distribution of random coefficients
f = zeros(N,m); z = f;
for i=1:m
    [F, X] = ecdf(V(:,i)); [X,ind] = unique(X); F = F(ind);
    Z = (X - mean(X))/std(X); 
    z(:,i) = linspace(min(Z),max(Z),N)';
    f(:,i) = pchip(Z,F,z(:,i));
end

%% Karhunen-Loeve Expansion
Prnd = repmat(Pmu,1,nrnd) + U*S(1:m,1:m)*randn(m,nrnd);
