% Autoencoder for 2-sphere
clc; close all; clearvars; 

m = 3; n = 3; l = 1; N = 1000;
% random ambient samples
C = 1 - 2*rand(N,m); Cnew = 1 - 2*rand(5000,m);
% sphere data
% X = repmat(1./sqrt(sum(C.^2,2)),1,m).*C; Xnew = repmat(1./sqrt(sum(Cnew.^2,2)),1,m).*Cnew;
% cube data
X = C./repmat(max(abs(C),[],2),1,m); Xnew = Cnew./repmat(max(abs(Cnew),[],2),1,m);

fig = figure;
scatter3(X(:,1), X(:,2), X(:,3),'o')
axis equal; fig.CurrentAxes.Visible = 'off';
% initialize with random entries
W = rand(n, m); c= rand(n,1); b = rand(m,1);
% normalize W
W = W./sqrt(sum(W.^2, 2));
L0 = loss(X,W,c,b);

rng default; % For reproducibility
gs = GlobalSearch('Display','iter');
problem = createOptimProblem('fmincon','x0',[reshape(W,1,m*n), c', b'],...
          'objective',@(theta) loss(X,reshape(theta(1:m*n),n,m),...
                      theta(m*n+1:m*n+n)',...
                      theta(m*n+n+1:end)'),...
          'lb',-100*ones(m*n+n+m,1),'ub',100*ones(m*n+n+m,1));
[theta_opt, loss_opt] = run(gs,problem);
Wopt = reshape(theta_opt(1:m*n),n,m);
copt = theta_opt(m*n+1:m*n+n)';
bopt = theta_opt(m*n+n+1:end)';

[Lopt,Xtilde,Hopt] = loss(Xnew,Wopt,copt,bopt);


hold on;
scatter3(Xtilde(:,1),Xtilde(:,2),Xtilde(:,3),'.')

function [L,Xtilde,H] = loss(X,W,c,b)
    % basis
%     sig = @(z) 1./(1 + exp(-z));
    sig = @(z) z;
    % encoder
    H = sig(X*W' + repmat(c', size(X,1), 1));
    % decoder
    Xtilde = sig(H*W + repmat(b', size(X,1), 1));
    L = sum(sqrt(sum((X - Xtilde).^2,2)));
end
