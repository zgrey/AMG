% Active Subspace approximation demo...
% Zach Grey, 02/12/2020
%
% The following implementations are based on the code and work developed by
% Paul Constantine and his research group including but not limited to:
% Zach Grey
% Andrew Glaws
% Jeffery Hokanson
% Izabel Aguiar
% Kerrek Stinson
% Paul Diaz
%
% Instructions: assuming you have pairs of inputs/outputs (x0,f(x0)) and a 
% uniform domain over a hypercube you'll want to specify the following:
% ub: a row vector of upper bounds
% lb: a row vector of lower bounds
% X0: a N by m matrix of inputs (in physical scales) 
% F:  a N by 1 column vector of outputs
% m:  the total number of parameters being varried (can be inferred from
%     the size of X0)
% N:  the total number of samples from the domain (can be inferred from the
%     size of X0)
clc; close all; rng(42);
% assign this to your local directory containing the matlab scripts
% codes available at GitHub, run the following:
% >>git clone https://github.com/zgrey/active_subspaces.git
% Linux
% AS_HOME = '/local/tmp/active_subspaces/matlab/';
% Windows
AS_HOME = 'C:\Users\zgrey\Documents\GitHub\active_subspaces\matlab\';

% add routines for training monomials
addpath([AS_HOME,'ResponseSurfaces'])
addpath([AS_HOME,'Domains'])

%% Dummay data from uniform tensor product domain
%%%%%%%%%%%%%%%%%%%%%%%%%%% THINGS TO MODIFY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total number of parameters (m ~= 24)
m = 10;
% number of random samples
N = 1000;
% given upper and lower bounds (replace these with your own domain def.):
ub = 2*ones(1,m); lb = ones(1,m);
% generate uniform random samples from box domain (replace with your data)
X0 = repmat(lb,N,1) + rand(N,m)*diag(ub - lb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define an affine transformation x = M*x0 + b so that x in [-1,1]^m
M = diag(2./(ub - lb)); b = -M*mean([ub; lb],1)';
% the scaled samples become
X = X0*M + repmat(b',N,1);
% and transforming back becomes (these are the same as X0)
X00 = (X - repmat(b',N,1))*diag(diag(1./M));

%% Dummy functions representing a response
% some random ridge directions
A = 2*rand(m,2) -1; A = A*diag(1./sqrt(sum(A.^2,1)));
% A(:,1) = [1; zeros(m-1,1)]; A(:,2) = [zeros(m-1,1); 1];
% a quadratic ridge function
f = @(X0) X0*A(:,1) + sum(X0*A(:,2)*A(:,2)'.*X0,2);
% a nonlinear (transcendental) ridge function
% f = @(X0) cos(2*pi*X0*A(:,1)) + sin(pi*X0*A(:,1)).^2 + (X0*A(:,1)) + (X0*A(:,1).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% THINGS TO MODIFY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate the function (replace this with your vector of responses)
F = f(X0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global linear model (strictly reduces to 1-dimensional subspace)
% fit a global linear model (make sure to use the transformed parameters!)
u = [ones(N,1), X] \ F;
% make this an orthonormal basis
u = u(2:end)/norm(u(2:end));
% make the shadow plot
scatter(X*u,F,50,'filled'); alpha(0.5); hold on;

% determine the extent of the original domain projected to the new 1d-domain
ext1 = linprog(u,[],[],[],[],-ones(m,1),ones(m,1));
ext2 = linprog(-u,[],[],[],[],-ones(m,1),ones(m,1));
% train a polynomial surrogate
p = 5; [Coef, ~, ~, ~, f_hat, res] = poly_train(X*u, F, p);
Rsqd = 1 - cov(res)./cov(F);
% reevaluate and plot the low-dimensional polynomial
Ydata = X*u; Yd = linspace(min(Ydata),max(Ydata),100)';
Yext = [ext1'; ext2']*u; Ye = linspace(min(Yext),max(Yext),100)';
[Pext,dPext] = poly_predict(Ye,Coef,p);
[P,dP] = poly_predict(Yd,Coef,p);
% plot the resulting approximation
plot(Yd,P,'k','linewidth',2); ax = gca; Ylim = ax.YLim; 
plot(Ye,Pext,'k--','linewidth',2); ax.YLim = Ylim;
% plot the extent of the original domain
scatter([ext1'; ext2']*u,min(ax.YLim)*ones(1,2),50,'filled','k');
% derivatives (if you're interested)
% quiver(Y(1:end-1),P(1:end-1),diff(Y),dP(1:end-1).*diff(Y),1); axis equal;
title(['Shadow Plot - Global Linear Model, R^2 = ', num2str(Rsqd)]);

%% Global quadratic model (reduces to arbitrary n-dimensional subspace)
% fit global quadratic model (make sure to use the transformed parameters!)
if N < nchoosek(m + 2,2)
    disp('ERROR: Not enough evaluations for quadratic model');
    fprintf('You need at least %f more samples', nchoosek(m + 2,2) - N + 1)
    fprintf('We recommend between %f and %f total samples as a rule of thumb...'...
            ,2*nchoosek(m + 2,2),10*nchoosek(m + 2,2))
else
    [~, ~, a, H, ~, ~] = poly_train(X,F,2);
    [W, Lambda] = eig(H*H + a*a'); Lambda = abs(Lambda);
    [Lambda, ind] = sort(diag(Lambda), 1, 'descend');
    W = W(:, ind);
end

% make shadow plots
figure;
% make a 1-dimensional shadow plot
subplot(1,2,1), scatter(X*W(:,1),F,50,'filled'); alpha(0.5); hold on;
% make a 1-dimensional surrogate
% determine the extent of the original domain projected to the new 1d-domain
ext1 = linprog(W(:,1),[],[],[],[],-ones(m,1),ones(m,1));
ext2 = linprog(-W(:,1),[],[],[],[],-ones(m,1),ones(m,1));
% train a polynomial surrogate
p = 3; [Coef, ~, ~, ~, f_hat, res] = poly_train(X*W(:,1), F, p);
Rsqd = 1 - cov(res)./cov(F);
% reevaluate and plot the low-dimensional polynomial
Ydata = X*u; Yd = linspace(min(Ydata),max(Ydata),100)';
Yext = [ext1'; ext2']*u; Ye = linspace(min(Yext),max(Yext),100)';
[Pext,dPext] = poly_predict(Ye,Coef,p);
[P,dP] = poly_predict(Yd,Coef,p);
% plot the resulting approximation
plot(Yd,P,'k','linewidth',2); ax = gca; Ylim = ax.YLim; 
plot(Ye,Pext,'k--','linewidth',2); ax.YLim = Ylim;
% plot the extent of the original domain
scatter([ext1'; ext2']*u,min(ax.YLim)*ones(1,2),50,'filled','k');
title(['1d Shadow Plot - Global Quadratic Model, R^2 = ', num2str(Rsqd)]);

% make a 2-dimensional shadow plot
subplot(1,2,2), scatter(X*W(:,1),X*W(:,2),50,'filled','cdata',F);
alpha(0.8); hold on; axis equal;
% make a 2-dimensional surrogate
% compute vertices of the zonotope
[~,ZV] = zonotope_vertices(W(:,1:2),10); th = atan2(ZV(:,2),ZV(:,1));
ZV = sortrows([th,ZV]); ZV = ZV(:,2:end);
Ydata2 = X*W(:,1:2); [Y1d2,Y2d2] = meshgrid(linspace(min(ZV(:,1)),max(ZV(:,1)),100)',...
                                  linspace(min(ZV(:,2)),max(ZV(:,2)),100)');
YZ = [reshape(Y1d2,100^2,1), reshape(Y2d2,100^2,1)];
X2 = YZ*W(:,1:2)';
% determine points within full extent of domain                                
indOUT = (X2 < -1 | X2 > 1); ind = max(indOUT,[],2) == 1;
% determine points within zonotope
ind = in2DZ(ZV,YZ);

[Coef2d, ~, ~, ~, f_hat2d, res2d] = poly_train(X*W(:,1:2), F, p);
Rsqd2d = 1 - cov(res2d)./cov(F);
[P2,dP2] = poly_predict([reshape(Y1d2,100^2,1), reshape(Y2d2,100^2,1)],Coef2d,p);

scatter(ZV(:,1),ZV(:,2),50,'k','filled'); plot([ZV(:,1); ZV(1,1)],[ZV(:,2); ZV(1,2)],'k','linewidth',2);
P2(~ind) = NaN;
contour(Y1d2,Y2d2,reshape(P2,100,100),50,'linewidth',1)
title(['2d Shadow Plot - Global Quadratic Model, R^2 = ' num2str(Rsqd2d)]); colorbar;

%% Evaluate subspace predictions
figure;
% plot eigenvalues
stem(Lambda,'filled'); grid on; ax = gca; ax.YScale = 'log';
title 'Eigenvalues'
figure;
% plot eigenvectors
stem(W(:,1:2),'filled'); grid on; hold on;
title Eigenvectors; legend('1st eig. vec.','2nd eig. vec.')