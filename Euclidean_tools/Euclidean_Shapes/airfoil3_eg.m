% turbine airfoil 3D embedding example
clc; close all; clearvars;
load turbine_airfoil.mat;
NN = 250;

% 4-point circle test case
PP = [cos(0:0.25*2*pi:2*pi)',sin(0:0.25*2*pi:2*pi)'];

%% 3D circular embedding
[PP3,n3,curv,a_spl,s_spl,s,alpha,temb,PP0,rad_ind,ttemb] = embrep3(PP,NN,'uni'); 

%% Visualize shape
set(0,'defaulttextInterpreter','latex','DefaultAxesFontSize',16)
fig = figure; hold on; axis equal;
% original points
% quiver(zeros(size(PP,1),1),zeros(size(PP,1),1),PP0(1:end-1,1),PP0(1:end-1,2),1,'color',0.75*[1 1 1],'maxheadsize',0.1);
% plot(PP(:,1),PP(:,2),'ko','markersize',5);
% 3D circular embedding
h = plot(PP3(:,1),PP3(:,2),'o-','linewidth',1);
% interpolated original points
% plot(PP0(:,1),PP0(:,2),'k.');
% unit normals
quiver(PP3(:,1),PP3(:,2),n3(:,1),n3(:,2),1,'color',h.Color);
axis equal; fig.CurrentAxes.Visible = 'off';
% curvature
fig = figure; hold on; axis equal;
bndplot(fig,PP3(:,1),PP3(:,2),log(curv)); fig.CurrentAxes.Visible = 'off';
plot(PP(:,1),PP(:,2),'ro','markersize',5);
% plot(PP3(:,1),PP3(:,2),'k.');
% caxis([-6 4]);
% alpha(t)
figure; plot(s-pi,alpha,'ko','markersize',4); grid on; hold on;
plot(ppval(s_spl,linspace(0,1,NN))-pi,ppval(a_spl,linspace(0,1,NN)),'k','linewidth',1);
xlabel('$$\theta$$','fontsize',20); ylabel('$$\hat{\alpha}(\theta)$$','fontsize',20);
% plot component functions alpha(t) and s(t)
figure;
plot(linspace(0,1,NN),ppval(a_spl,linspace(0,1,NN))); hold on;
plot(linspace(0,1,NN),1/pi*ppval(s_spl,linspace(0,1,NN)));

%% regular component fucntion spline
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
legend('\alpha','s','x','y')
% curvature of regular spline over airfoil
dx = 1/(ub-lb)*ppval(fnder(x_spl,1),tt); ddx = 1/(ub-lb)^2*ppval(fnder(x_spl,2),tt);
dy = 1/(ub-lb)*ppval(fnder(y_spl,1),tt); ddy = 1/(ub-lb)^2*ppval(fnder(y_spl,2),tt);
curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);
fig = figure; axis equal;
bndplot(fig,ppval(x_spl,tt),ppval(y_spl,tt),log(curv2)); fig.CurrentAxes.Visible = 'off';
plot(PP(:,1),PP(:,2),'ro','markersize',5);
% caxis([-6 4]);
% compare curvatures
figure;
semilogy(ttemb,curv,'linewidth',2); hold on; grid on;
semilogy(tt,curv2,'--','linewidth',2,'color',0.75*ones(1,3));
legend('emb. rep.','reg. spline');