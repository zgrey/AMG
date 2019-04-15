% ralphie 3D embedding example
clc; close all; clearvars;
skp = 10; NN = 250; G = [0 -1;1 0];

%% image processing
Rim = imread('./ralphie.png'); BW = imbinarize(Rim(:,:,1));
Rbnd = bwboundaries(~BW); bnd = Rbnd{3}*G; % use Rbnd{2}-inner or Rbnd{3}-outer
% select subset (ensure point set is from boundary of a simply-connected domain)
bndsub = bnd(1:skp:end,:); N = size(bndsub,1);
% shift and scale
bndsub = 2*(bndsub - min(bndsub,[],1))./(max(bndsub,[],1)-min(bndsub,[],1)) - 1;
x0 = [mean(bndsub(:,1)),mean(bndsub(:,2))]; bndsub = bndsub - x0;

%% 3D circular embedding
[PP3,n3,curv,a_spl,s_spl,s,alpha,t] = embrep3(bndsub,NN,'curv'); 

%% Visualize shape
fig = figure; hold on; axis equal;
plot([bndsub(:,1);bndsub(1,1)],[bndsub(:,2);bndsub(1,2)],'k','linewidth',5,'Color',[0.929,0.694,0.125]);
% 3D circular embedding
plot(PP3(:,1),PP3(:,2),'k','linewidth',2);
quiver(PP3(:,1),PP3(:,2),-n3(:,1),-n3(:,2),1,'k');
scatter(bndsub(:,1),bndsub(:,2),'k');
axis equal; fig.CurrentAxes.Visible = 'off';
% curvature
fig = figure; hold on; axis equal;
bndplot(fig,PP3(:,1),PP3(:,2),log(curv)); fig.CurrentAxes.Visible = 'off';
% 3D embedding
figure; plot3(t,s,alpha,'o'); grid on; hold on;
plot3(linspace(0,1,NN),ppval(s_spl,linspace(0,1,NN)),ppval(a_spl,linspace(0,1,NN)),'k','linewidth',2);
xlabel '$$t$$'; ylabel '$$\hat{\theta}(t)$$'; zlabel '$$\hat{\alpha}(t)$$';
% plot component functions alpha(t) and s(t)
figure;
plot(linspace(0,1,NN),ppval(a_spl,linspace(0,1,NN))); hold on;
plot(linspace(0,1,NN),1/pi*ppval(s_spl,linspace(0,1,NN)));

%% regular spline comparison
PP = bndsub;
% repeat first point to close shape
P = [PP; PP(1,:)];
% discrete length scales for splines
t0 = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 + ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
ub = max(t0); lb = 0;
% scale to [0,1]
t = 1/(ub-lb)*t0;
x_spl = csape(t,P(:,1),'periodic');
y_spl = csape(t,P(:,2),'periodic');
% plot component functions
tt = linspace(0,1,NN)';
plot(tt,ppval(x_spl,tt)); hold on;
plot(tt,ppval(y_spl,tt));
legend('$$\alpha$$','s','x','y')
% curvature of regular spline
dx = 1/(ub-lb)*ppval(fnder(x_spl,1),tt); ddx = 1/(ub-lb)^2*ppval(fnder(x_spl,2),tt);
dy = 1/(ub-lb)*ppval(fnder(y_spl,1),tt); ddy = 1/(ub-lb)^2*ppval(fnder(y_spl,2),tt);
curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);
fig = figure; axis equal;
bndplot(fig,ppval(x_spl,tt),ppval(y_spl,tt),log(curv2)); fig.CurrentAxes.Visible = 'off';
figure;
semilogy(tt,curv); hold on; grid on;
semilogy(tt,curv2,'--');
legend('emb. rep.','reg. spline');