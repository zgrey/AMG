% Implicit RBF Surface

function [P,Phi,GP,G] = rbf_surf(U,w,M,S)
% U in R^(N x d) is a stacked collection of R^d vectors of the surf domain
% w in R^(m) are the basis coefficients
% M in R^(m x d) are the centers of the RBF's
% S in R^(d x d) defines the uniform length scales and eccentricity

m = size(M,1); N = size(U,1); d = size(U,2);

Phi = zeros(N,m); G   = zeros(N*d,m);
parfor j=1:m
    % Gaussian Kernel
    D = bsxfun(@minus,U,M(j,:));
    Phi(:,j)  = exp(-sum(D.*(D*S),2));
    
    % Gradient 
    G(:,j)    = reshape(bsxfun(@times,Phi(:,j),-2*D*S),N*d,1);
end

P  = [U, Phi*w]; GP = reshape(G*w,N,d);

% Note, this is how to reshape data for Matlab contours/surfaces/meshes
% X = reshape(P(:,1),sqrt(N),sqrt(N));
% Y = reshape(P(:,2),sqrt(N),sqrt(N));
% Z = reshape(P(:,3),sqrt(N),sqrt(N));