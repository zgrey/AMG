function [W, lambda] = fasshauer_periodic_kernel(x, n, beta)
%FASSHAUER_PERIODIC_KERNEL
%
% Return eigenpairs of periodic kernel from Appendix A.5.2 of Fasshauer and
% McCourt, 'Kernel-based Approximation Methods using MATLAB'
%
% INPUTS
% x         grid 
% n         number of eigenpairs
% beta      controls smoothness, larger beta is smoother
%
% OUTPUTS
% W         evaluations of first n eigenvecs at grid points in x
% lambda    first n eigenvals
%
% Example:
% % Set up a grid covering one period
% nx = 101; x = linspace(0, 1, nx)';
% 
% % Compute eigenpairs
% ne = 100; beta = 2;
% [W, lam] = fasshauer_periodic_kernel(x, ne, beta);
%
% % Draw random Gaussians
% nsamples = 10; z = randn(ne, nsamples);
%
% % Evaluate the zero-mean Karhunen-Loeve expansion for each sample
% F = W*(bsxfun(@times, sqrt(lam), z));
%

% make column vec
x = x(:);

% eigenvals
j = (1:ceil(n/2));
lambda = [(2*pi*j).^(-2*beta); (2*pi*j).^(-2*beta)];
lambda = lambda(:); lambda = lambda(1:n);

% eigenvecs
W = zeros(length(x), 2*ceil(n/2));
W0 = sqrt(2)*sin(2*pi*x*j);
W1 = sqrt(2)*cos(2*pi*x*j);
W(:,1:2:end) = W0; W(:,2:2:end) = W1;
W = W(:,1:n);

