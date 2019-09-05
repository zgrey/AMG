% ralphie 3D embedding example
clc; close all; clearvars;

%% Number of refined samples
NN = 25;

%% image processing
skp = 4; Rim = imread('./ralphie.png');
% skp = 10; Rim = imread('./MPEG7dataset/original/octopus-16.gif');
% skp = 10; Rim = imread('./MPEG7dataset/original/bat-1.gif');
% skp = 20; Rim = imread('./MPEG7dataset/original/dog-1.gif');
% skp = 10; Rim = imread('./MPEG7dataset/original/beetle-6.gif');

BW = imbinarize(Rim(:,:,1)); Rbnd = bwboundaries(~BW);

bnd = Rbnd{2}*[0 -1;1 0];

% flip to be counter-clockwise
bnd = bnd(end:-1:1,:);

% select subset (ensure point set is from boundary of a simply-connected domain)
bndsub = bnd(1:skp:end,:); N = size(bndsub,1);
% shift and scale
bndsub = 2*(bndsub - min(bndsub,[],1))./(max(bndsub,[],1)-min(bndsub,[],1)) - 1;
x0 = [mean(bndsub(:,1)),mean(bndsub(:,2))]; bndsub = bndsub - x0;

%% 3D circular embedding
% get default options
options = embOptions([]); 
% options.AffineTrans = 'none';
% compute embedding representation
emb = embrep3(bndsub,NN,'curv',options);

%% Visualize shape
set(0,'defaulttextInterpreter','latex','DefaultAxesFontSize',16)
fig = figure; hold on; axis equal;
% plot nominal landmarks
plot(bndsub(:,1),bndsub(:,2),'ko','markersize',5);
% plot embedding representation
plot(emb.pts(:,1),emb.pts(:,2),'linewidth',1);
% plot unit normals
figure(1);
quiver(emb.pts(:,1),emb.pts(:,2),emb.nml(:,1),emb.nml(:,2),1);
axis equal; fig.CurrentAxes.Visible = 'off';

% curvature
fig = figure; hold on; axis equal;
bndplot(fig,emb.pts(:,1),emb.pts(:,2),log(emb.curv)); fig.CurrentAxes.Visible = 'off';
plot(emb.pts(:,1),emb.pts(:,2),'k.');
caxis([-6 4]);

% plot inner products at nominal length scales
figure; plot(emb.TF.t,emb.TF.alph,'ko','markersize',5); grid on; hold on;
% plot inner products at refined length scales
plot(emb.t,emb.alph.eval,'k','linewidth',1);
xlabel('$$t$$','fontsize',20); ylabel('$$\hat{\alpha}(t)$$','fontsize',20);

% plot shifted angles (for visualization) at nominal length scales
figure; plot(emb.TF.t,emb.TF.th,'kd','markersize',5); grid on; hold on;
% plot shifted angles (for visualization) at refined length scales
plot(emb.t,emb.th.eval,'k','linewidth',1);
xlabel('$$t$$','fontsize',20); ylabel('$$\hat{\theta}(t)$$','fontsize',20);

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
tt = linspace(0,1,NN)';
% recompute points
P = [ppval(x_spl,tt), ppval(y_spl,tt)];
% recompute discrete lengths
L = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 +...
            ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
        
% curvature of regular spline
dx = 1/(ub-lb)*ppval(fnder(x_spl,1),tt); ddx = 1/(ub-lb)^2*ppval(fnder(x_spl,2),tt);
dy = 1/(ub-lb)*ppval(fnder(y_spl,1),tt); ddy = 1/(ub-lb)^2*ppval(fnder(y_spl,2),tt);
curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);
% curvature at nominal landmarks
dx = 1/(ub-lb)*ppval(fnder(x_spl,1),t); ddx = 1/(ub-lb)^2*ppval(fnder(x_spl,2),t);
dy = 1/(ub-lb)*ppval(fnder(y_spl,1),t); ddy = 1/(ub-lb)^2*ppval(fnder(y_spl,2),t);
nom_curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);

% visualize curvature
fig = figure; axis equal;
bndplot(fig,ppval(x_spl,tt),ppval(y_spl,tt),log(curv2)); fig.CurrentAxes.Visible = 'off';
plot(ppval(x_spl,tt),ppval(y_spl,tt),'k.');
caxis([-6 4]);
% compare curvatures
figure;
h1 = semilogy(emb.t,emb.curv); hold on; grid on;
h2 = semilogy(tt,curv2,'--');
semilogy(emb.nom.t,emb.nom.curv,'*','color',h1.Color);
semilogy(t,nom_curv2,'o','color',h2.Color)
legend('emb. rep.','reg. spline');

%% build radially convex approximation
% compute monotonic indices (remove last repeated point)
if emb.nom.th(1) < emb.nom.th(end)
    incr = mono_vec(emb.nom.th(1:end-1),1e-3);
else
    incr = mono_vec(-emb.nom.th(1:end-1),1e-3);
end
% 3D circular embedding
% use consistent affine transformation
options.AffineTrans = 'custom'; 
options.M = emb.TF.M; options.b = emb.TF.b; options.Minv = emb.TF.Minv;
% use monotonic interpolation
options.ThetaSpline = 'pchip';
star_emb = embrep3(emb.nom.pts(incr,:),NN,'uni',options);

%% Visualize shape
set(0,'defaulttextInterpreter','latex','DefaultAxesFontSize',16)
fig = figure; hold on; axis equal;
% plot radially convex embedding representation
plot(star_emb.pts(:,1),star_emb.pts(:,2),'linewidth',2);
% plot nominal landmarks
plot(bndsub(:,1),bndsub(:,2),'ko','markersize',5);
% plot unit normals
% quiver(star_emb.pts(:,1),star_emb.pts(:,2),star_emb.nml(:,1),star_emb.nml(:,2),1,'k');
axis equal; fig.CurrentAxes.Visible = 'off';

% curvature
fig = figure; hold on; axis equal;
bndplot(fig,star_emb.pts(:,1),star_emb.pts(:,2),log(star_emb.curv)); fig.CurrentAxes.Visible = 'off';
plot(star_emb.pts(:,1),star_emb.pts(:,2),'k.');
caxis([-6 4]);

% plot inner products at nominal length scales
figure; plot(star_emb.nom.t,star_emb.nom.alph,'ko','markersize',5); grid on; hold on;
% plot inner products at refined length scales
plot(star_emb.t,star_emb.alph.eval,'k','linewidth',1);
xlabel('$$t$$','fontsize',20); ylabel('$$\hat{\alpha}(t)$$','fontsize',20);

% plot shifted angles (for visualization) at nominal length scales
figure; plot(star_emb.nom.t,star_emb.nom.th + pi,'kd','markersize',3); grid on; hold on;
% plot shifted angles (for visualization) at refined length scales
plot(star_emb.t,star_emb.th.eval + pi,'k','linewidth',1);
xlabel('$$t$$','fontsize',20); ylabel('$$\hat{\theta}(t)$$','fontsize',20);