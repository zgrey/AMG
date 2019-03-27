% Riemannian view on shape optimization [2014]
% V. Schulz numerical example
close all; clear variables; clc;

N     = 500;  % Number of discrete shape "particles"
mu    = 2;    % function parameter (see [Schulz, 2014])
alpha = 0.75; % crude, fixed relaxation factor (not exact line search, per paper)

% Schulz's initial shape
s = linspace(0,2*pi,N)';
c0 = 1/2*[cos(s)-0.15*abs(1-sin(2*s)).*cos(s),...
          sin(s)-0.15*abs(1-cos(2*s)).*cos(s)];

% Random initial shape (WIP, requires rbf_surf.m)
% addpath ../SHDP/
% NN = 10;
% e  = 0.05;
% pp = 0;
% slpe  = 0.5;
% c0 = rndshapes(pp,slpe,e,NN,N);

% plot nominal shape
plot(c0(:,1),c0(:,2)); axis equal; hold on; title 'Initial Shape'

% Discrete normal calculation (central diff.)
G = [0 1;-1 0];
n = (c0(2:end,:)-c0(1:end-1,:))*G';
n = n./sqrt(sum(n.^2,2));

% plot surface normals
cc = (c0(2:end,:)+c0(1:end-1,:))/2;
quiver(cc(:,1),cc(:,2),n(:,1),n(:,2));

% shape distance approximation
d = @(c) 1/(2*mu)*sum(atan2(c(2:end,2),c(2:end,1))-atan2(c(1:end-1,2),c(1:end-1,1)).*...
    (abs(sqrt(c(2:end,1).^2 + c(2:end,2).^2)-1) + ...
    abs(sqrt(c(1:end-1,1).^2 + c(1:end-1,2).^2)-1)));

% Newton's method
figure; plot(c0(:,1),c0(:,2)); axis equal; hold on; title Iterates
c = cc; gradf = ones(N,1); i = 1;
while norm(gradf) > 1e-7
    gradf = c(:,1).^2 + mu^2*c(:,2).^2-1;
    Hessf = 2*sqrt(c(:,1).^2 + mu^4*c(:,2).^2);
    R = -1./Hessf.*gradf;
    c = c + alpha*R.*n;
    
    CI(i) = norm(gradf); i = i + 1;
    
    plot([c(:,1);c(1,1)],[c(:,2);c(1,2)]);
end

% more plots and linear convergence rate calculation
figure; plot([c(:,1);c(1,1)],[c(:,2);c(1,2)],'k','linewidth',3);
axis equal; title 'Final Shape'
figure; semilogy(1:i-1,CI,'k-o','linewidth',2); grid on; title Convergence;
ylabel '||grad(f)||'
A = [ones(i-1,1) log([1:i-1])']; u = A'*A \ A'*CI';

% Results from paper
figure; semilogy([1:5],[0.9222,0.1382,0.3571e-2,0.8187e-5,0.1736e-9],'k-o','linewidth',2);
grid on; title 'Reported Convergence - Table 1 [Schulz, 2014]'; ylabel 'd(c)';