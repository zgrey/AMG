% an example using self-intersection detection routine
clc; clearvars; close all;

n = 500;
% parametrize curve
t = linspace(-1,1,n)';
% Leminscate with one known intersection
X = [cos(pi*t), sin(pi*t).*cos(pi*t)]./(1+sin(pi*t).^2);

% run routine
tic;
[~,Xint,ind] = selfintersect(X);
runtime = toc;
disp('Runtime (sec):')
disp(runtime);
disp('Total intersections found:')
disp(sum(sum(ind)));
tic; [~,p,q] = step_selfintersect(X); runtime = toc; ind = [p,q];
disp('Step runtime (sec):')
disp(runtime);

% plot results
plot(X(:,1), X(:,2),'LineWidth',2); hold on; axis equal;
scatter(Xint(:,1), Xint(:,2),100,'ro','linewidth', 2);
scatter(X(ind,1),X(ind,2),25,'ro','filled')
quiver(X(ind,1), X(ind,2), X(ind+1,1) - X(ind,1), X(ind+1,2) - X(ind,2), 0,'r','LineWidth',2);

%% try with some crazy CST shapes
addpath ../../../ASAP/AUTO/;
% optimize to find nominal values
m = 100; nl = ceil(n/2); 
% an "intuitive" shape sampling (clustered at LE and TE)
l = (cos(linspace(0,pi,nl)) + 1)/2;

dv0 = 0.5*[ones(m/2,1); ones(m/2,1)];
% Define upper and lower bounds
pct = 1;
% lower surface
lb0(1:m/2) = (1+pct)*dv0(1:m/2); ub0(1:m/2) = (1-pct)*dv0(1:m/2);
% upper surface
ub0(m/2+1:m) = (1+pct)*dv0(m/2+1:end); 
lb0(m/2+1:m) = (1-pct)*dv0(m/2+1:end);
% random hypercube samples
X = 2*rand(1,length(dv0))-1;
dv = bsxfun(@plus,lb0,bsxfun(@times,ub0-lb0,0.5*(X(1,:)+1)));
[coordU, coordL] = cst_airfoil(l',dv(1:m/2)',dv(m/2+1:end)',0);
bnd = [coordL'; coordU(:,end:-1:1)'];
[bnd,Ntru] = unique_points(bnd); bnd = [bnd; bnd(1,:)];

% run routine
tic;
[~,Xint,ind] = selfintersect(bnd);
runtime = toc;
disp('Runtime (sec):')
disp(runtime);
disp('Total intersections found:')
disp(sum(sum(ind)));
tic; [~,p,q] = step_selfintersect(bnd); runtime = toc; ind = [p,q];
disp('Step runtime (sec):')
disp(runtime);

% plot results
figure;
plot(bnd(:,1), bnd(:,2),'LineWidth',1); hold on; axis equal;
scatter(Xint(:,1), Xint(:,2),100,'ro','linewidth', 2);
scatter(bnd(ind,1),bnd(ind,2),25,'ro','filled')
quiver(bnd(ind,1), bnd(ind,2),bnd(ind+1,1) - bnd(ind,1), bnd(ind+1,2) - bnd(ind,2),1,'r','LineWidth',2);

%% try using the graph of a polynomial with known roots
Nrts = 5; npts = 500;
rts = linspace(0.1,0.9,Nrts); x = linspace(0,1,npts)';
coef = -10 + 10*rand(Nrts,1);
poly = ones(npts,1);
for i=1:Nrts
    poly = coef(i)*poly.*(x - rts(i));
end
grph = [[0,0]; [x, poly]; [1,0]; [0,0]];

% run routine
tic;
[~,Xint,ind] = selfintersect(grph);
runtime = toc;
disp('Runtime (sec):')
disp(runtime);
disp('Total intersections found:')
disp(sum(sum(ind)));
tic; [~,p,q] = step_selfintersect(grph); runtime = toc; ind = [p,q];
disp('Step runtime (sec):')
disp(runtime);

% plot results
figure;
plot(grph(:,1), grph(:,2), 'LineWidth', 2); hold on;
scatter(Xint(:,1), Xint(:,2),100,'ro','linewidth', 2);
scatter(grph(ind,1),grph(ind,2),25,'ro','filled')
disp('Worst-case error:')
disp(max(abs(Xint(:,1) - rts')))

%% try using the graph of sine
freq = 15; npts = 500;
x = linspace(0,1,npts)';
sine = sin(pi*x*freq);
grph = [[0,0]; [x, sine]; [1,0]; [0,0]];

% run routine
tic;
[bool,Xint,ind] = selfintersect(grph);
runtime = toc;
disp('Runtime (sec):')
disp(runtime);
disp('Total intersections found:')
disp(sum(sum(ind)));
tic; [~,p,q] = step_selfintersect(grph); runtime = toc; ind = [p,q];
disp('Step runtime (sec):')
disp(runtime);

% plot results
figure;
plot(grph(:,1), grph(:,2), 'LineWidth', 2); hold on;
scatter(Xint(:,1), Xint(:,2),100,'ro','linewidth', 2);
scatter(grph(ind,1),grph(ind,2),25,'ro','filled')