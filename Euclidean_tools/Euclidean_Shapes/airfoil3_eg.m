% turbine airfoil 3D embedding example
clc; close all; clearvars; rng(42);
load turbine_airfoil.mat;
NN = 250;

% 4-point circle test case (coordinate aligned)
% PP = [cos(0:0.5*pi:2*pi)',sin(0:0.5*pi:2*pi)'];
% tru_crv = ones(NN,1);

% 4-points circle test case (coordinate misaligned)
% G = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
% PP = [cos(0:0.5*pi:2*pi)',sin(0:0.5*pi:2*pi)']*G';
% tru_crv = ones(NN,1);

% peanut test case
% G = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
% scl = 1; Pshp =@(t) [scl*(cos(4*pi*t)+2).*cos(2*pi*t), (cos(4*pi*t)+2).*sin(2*pi*t)]*G';
% shp_curv = @(t) 24*sqrt(2)*scl*abs(-7-8*cos(4*pi*t)+cos(8*pi*t))./...
%                ( 18 + 34*scl^2 + (27+5*scl^2)*cos(4*pi*t) - 6*(5*scl^2-3)*cos(8*pi*t) + 9*cos(12*pi*t) - 9*scl^2*cos(12*pi*t) ).^(3/2);
% PP = Pshp(linspace(0,1,4)'); tru_crv = shp_curv(linspace(0,1,NN)');

% 4-point ellipse test case
% G = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
% scl = 2; PP = [scl*cos(0:0.25*2*pi:2*pi)',sin(0:0.25*2*pi:2*pi)']*G'; 
% tt = linspace(0,1,NN); tru_crv = scl./( scl^2*sin(2*pi*tt).^2 + cos(2*pi*tt).^2).^(3/2);

%% 3D circular embedding
% get default options
options = embOptions([]); 
% options.AffineTrans = 'none';
% compute embedding representation
emb = embrep3(PP,NN,'uni',options);
emb_interp_err = sum(sqrt(sum((PP - emb.nom.pts).^2,2)));
disp('Embedding Interpolation Verification:')
disp(['Interpolation error = ',num2str(emb_interp_err)])

%% Visualize shape
set(0,'defaulttextInterpreter','latex','DefaultAxesFontSize',16)
fig = figure; hold on; axis equal;
% plot nominal landmarks
plot(PP(:,1),PP(:,2),'ko','markersize',8);
% plot embedding representation
plot(emb.pts(:,1),emb.pts(:,2),'k','linewidth',1);
plot(emb.nom.pts(:,1),emb.nom.pts(:,2),'k*')
% plot unit normals
quiver(emb.pts(:,1),emb.pts(:,2),emb.nml(:,1),emb.nml(:,2),1,'k');
axis equal; fig.CurrentAxes.Visible = 'off';

% curvature
fig = figure; hold on; axis equal;
bndplot(fig,emb.pts(:,1),emb.pts(:,2),log(emb.curv)); fig.CurrentAxes.Visible = 'off';
plot(emb.pts(:,1),emb.pts(:,2),'w.','markersize',8);
% caxis([-6 4]);

% plot inner products at nominal length scales
figure; plot(emb.TF.t,emb.TF.alph,'ko','markersize',5); grid on; hold on;
% plot inner products at refined length scales
plot(emb.t,emb.alph.eval,'k','linewidth',1);
xlabel('$$t$$','fontsize',20); ylabel('$$\hat{\alpha}(t)$$','fontsize',20);

% plot shifted angles (for visualization) at nominal length scales
figure; plot(emb.TF.t,emb.TF.th + pi,'kd','markersize',3); grid on; hold on;
% plot shifted angles (for visualization) at refined length scales
plot(emb.t,emb.th.eval + pi,'k','linewidth',1);
xlabel('$$t$$','fontsize',20); ylabel('$$\hat{\theta}(t)$$','fontsize',20);

%% regular spline comparison
P = unique_points(PP);
% repeat first point to close shape
P = [P; P(1,:)];

% discrete length scales for splines
t0 = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 + ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
ub = max(t0); lb = 0;
% scale to [0,1]
t = 1/(ub-lb)*t0;
% t0 = 0:size(P,1)-1;
% t = t0;

% Build spline
% x_spl = csape(t,P(:,1),'periodic');
% y_spl = csape(t,P(:,2),'periodic');
xy_spl = csape(t,P','periodic');
tt = linspace(0,max(t),NN)';
% recompute points
% P = [ppval(x_spl,tt), ppval(y_spl,tt)];
P = fnval(xy_spl,tt)';
% recompute discrete lengths
L = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 +...
            ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
        
% curvature of regular spline
% dx = 1/(ub-lb)*ppval(fnder(x_spl,1),tt); ddx = 1/(ub-lb)^2*ppval(fnder(x_spl,2),tt);
% dy = 1/(ub-lb)*ppval(fnder(y_spl,1),tt); ddy = 1/(ub-lb)^2*ppval(fnder(y_spl,2),tt);
% curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);
d1 = fnval(fnder(xy_spl,1),tt)'; d2 = fnval(fnder(xy_spl,2),tt)';
dx = d1(:,1); ddx = d2(:,1);
dy = d1(:,2); ddy = d2(:,2);
curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);
% curvature at nominal landmarks
% dx = 1/(ub-lb)*ppval(fnder(x_spl,1),t); ddx = 1/(ub-lb)^2*ppval(fnder(x_spl,2),t);
% dy = 1/(ub-lb)*ppval(fnder(y_spl,1),t); ddy = 1/(ub-lb)^2*ppval(fnder(y_spl,2),t);
% nom_curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);
d1 = fnval(fnder(xy_spl,1),t)'; d2 = fnval(fnder(xy_spl,2),t)';
dx = d1(:,1); ddx = d2(:,1);
dy = d1(:,2); ddy = d2(:,2);
nom_curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);

tt0 = cumsum([0; sqrt(( PP(2:end,1) - PP(1:end-1,1) ).^2 + ( PP(2:end,2) - PP(1:end-1,2) ).^2)],1); tt0 = tt0/max(tt0);
reg_interp_err = sum(sqrt(sum((PP - fnval(xy_spl,tt0)').^2,2)));
disp('Regular Interpolation Verification:')
disp(['Interpolation error = ',num2str(reg_interp_err)])

%% visualize curvature
% regular spline
fig = figure; axis equal;
bndplot(fig,P(:,1),P(:,2),log(curv2)); fig.CurrentAxes.Visible = 'off';
plot(P(:,1),P(:,2),'w.','markersize',8);
% true curvature
% fig = figure; axis equal;
% if exist('tru_crv','var')
% bndplot(fig,PP(:,1),PP(:,2),log(tru_crv)); fig.CurrentAxes.Visible = 'off';
% plot(PP(:,1),PP(:,2),'w.','markersize',8);
% end
% caxis([-6 4]);
% compare curvatures
figure;
if exist('tru_crv','var')
semilogy(tt,tru_crv,'k','linewidth',2); hold on;
end
h2 = semilogy(tt,curv2,'--','linewidth',2);
h1 = semilogy(emb.t,emb.curv,'--','linewidth',2,'markersize',2); hold on; grid on;
semilogy(emb.nom.t,emb.nom.curv,'o','color',h1.Color,'linewidth',2);
semilogy(t,nom_curv2,'d','color',h2.Color,'linewidth',2)
xlabel '$$t$$'; ylabel '$$ \hat{\kappa}(t)$$';
% legend('emb. rep.','reg. spline');