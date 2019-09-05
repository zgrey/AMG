% convergence study
clc; close all; clearvars; rng(41);

%% airfoil study
addpath ../../../ASAP/AUTO/;
% optimize to find nominal values
m = 10; nl = 5001; 
% uniform over parameterization domain (not uniform over shape!)
% l = linspace(0,1,nl);
% an "intuitive" shape sampling (clustered at LE and TE)
l = (cos(linspace(0,pi,nl)) + 1)/2;
% an "unintuitive" shape sampling (clustered at mid-chord)
% l = (sin(linspace(pi/2,-pi/2,nl)) + 1)/2;
dv0 = 0.1*[-ones(m/2,1); ones(m/2,1)];
% Define upper and lower bounds
pct = 1;
% lower surface
lb0(1:m/2) = (1+pct)*dv0(1:m/2); ub0(1:m/2) = (1-pct)*dv0(1:m/2);
% upper surface
ub0(m/2+1:m) = (1+pct)*dv0(m/2+1:end); lb0(m/2+1:m) = (1-pct)*dv0(m/2+1:end);
% random hypercube samples
Nshp = 1;
X = 2*rand(Nshp,length(dv0))-1;

for ii = 1:Nshp
dv = bsxfun(@plus,lb0,bsxfun(@times,ub0-lb0,0.5*(X(ii,:)+1)));
[coordU, coordL] = cst_airfoil(l',dv(1:m/2)',dv(m/2+1:end)',0);
bnd = [coordL'; coordU(:,end:-1:1)'];
[bnd,Ntru] = unique_points(bnd); bnd = [bnd; bnd(1,:)];
tru_t = cumsum([0; sqrt(( bnd(2:end,1) - bnd(1:end-1,1) ).^2 + ( bnd(2:end,2) - bnd(1:end-1,2) ).^2)],1);
tru_t = tru_t/max(tru_t); 

%% shift and scale
x0 = [mean(bnd(:,1)),mean(bnd(:,2))]; bnd = bnd - x0;

Nmax = size(bnd,1);
skp = 2.^(0:floor(log(Nmax)/log(2))-2); skp = skp(end:-1:1);

%% Convergence study
for i=1:length(skp)
    % select subset (ensure point set is from boundary of a simply-connected domain)
    P0 = bnd(1:skp(i):end,:); 
    P0 = unique_points(P0); P0 = [P0; P0(1,:)]; Npts(i) = size(P0,1);
    P = P0;
    
    %% 3D circular embedding
    % get default options
    options = embOptions([]); 
%     options.AffineTrans = 'none';
    % compute embedding representation
    emb = embrep3(P0,250,'nom',options);
    
    %% regular spline comparison
    
    % discrete length scales for splines
    t0 = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 + ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
    ub = max(t0); lb = 0;
    
    % scale to [0,1]
    t = t0/ub;
    x_spl = csape(t,P(:,1),'periodic');
    y_spl = csape(t,P(:,2),'periodic');
    
    % recompute points at refined "true" landmarks
    P = [ppval(x_spl,tru_t), ppval(y_spl,tru_t)];
    % recompute discrete lengths
    L = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 +...
                ( P(2:end,2) - P(1:end-1,2) ).^2)],1);

    % curvature of regular spline
    dx = 1/(ub-lb)*ppval(fnder(x_spl,1),tru_t); ddx = 1/(ub-lb)^2*ppval(fnder(x_spl,2),tru_t);
    dy = 1/(ub-lb)*ppval(fnder(y_spl,1),tru_t); ddy = 1/(ub-lb)^2*ppval(fnder(y_spl,2),tru_t);
    curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);
    % curvature at nominal landmarks
    dx = 1/(ub-lb)*ppval(fnder(x_spl,1),t); ddx = 1/(ub-lb)^2*ppval(fnder(x_spl,2),t);
    dy = 1/(ub-lb)*ppval(fnder(y_spl,1),t); ddy = 1/(ub-lb)^2*ppval(fnder(y_spl,2),t);
    nom_curv2 = abs( dx.*ddy - dy.*ddx )./( dx.^2 + dy.^2 ).^(3/2);
 
    % reevaluate at points for error computation
    % convert between length scale measures
    treval = pchip(t,emb.TF.t,tru_t);
    emb_pts = (ppval(emb.alph.spl,treval).*[cos(ppval(emb.th.spl,treval)),sin(ppval(emb.th.spl,treval))] - repmat(emb.TF.b',length(tru_t),1))*emb.TF.Minv';
    reg_pts = P;
    % compute error
%     err_emb(i,ii) = 1/Nmax*sum(sqrt(sum((bnd - emb_pts).^2,2)));
%     err_reg(i,ii) = 1/Nmax*sum(sqrt(sum((bnd - reg_pts).^2,2)));
    
    % compute shape disctance over the Grassmannian
    bndTF = affine_trans(bnd(1:end-1,:),'LA');
    embTF = affine_trans(emb_pts(1:end-1,:),'LA');
    err_emb(i,ii) = dGr_np(1/sqrt(Ntru-1)*bndTF,1/sqrt(Ntru-1)*embTF);
    regTF = affine_trans(reg_pts(1:end-1,:),'LA');
    err_reg(i,ii) = dGr_np(1/sqrt(Ntru-1)*bndTF,1/sqrt(Ntru-1)*regTF);
    
%     plot(bndTF(:,1),bndTF(:,2),'k-o'); hold on;
%     plot(embTF(:,1),embTF(:,2),'g-*');
%     plot(regTF(:,1),regTF(:,2),'r-*');
    
    clc; disp(['AGGREGATING DATA...',num2str(i/length(skp)*100),'% Complete']);
    disp([num2str(ii),'/',num2str(Nshp),' Shapes Complete']);
end
end
%% convergence plots
set(0,'defaulttextInterpreter','latex','DefaultAxesFontSize',24)
figure; 

loglog(Npts,Npts.^(-4),'k','linewidth',2); hold on;

% Cartesian distances (errors)
h1 = loglog(Npts,mean(err_reg,2).^2,'--d','linewidth',2,'markersize',8);  grid on; hold on;

% embedding distances (errors)
h2 = loglog(Npts,mean(err_emb,2).^2,'-o','linewidth',2,'markersize',8);

err_max = max(err_reg,[],2).^2; err_min = min(err_reg,[],2).^2;
if min(err_min) == 0
    h = fill([Npts';Npts(end:-1:1)'],[err_max(1:end);1/2*eps;err_min(end-1:-1:1)],h1.Color);
    scatter(Npts(end),1/2*eps,'ko','linewidth',2)
else
    h = fill([Npts';Npts(end:-1:1)'],[err_max;err_min(end:-1:1)],h1.Color);
end
h.FaceAlpha = 0.15;
err_max = max(err_emb,[],2).^2; err_min = min(err_emb,[],2).^2;
if min(err_min) == 0
    h = fill([Npts';Npts(end:-1:1)'],[err_max(1:end);1/2*eps;err_min(end-1:-1:1)],h2.Color);
    scatter(Npts(end),1/2*eps,'ko','linewidth',2)
else
    h = fill([Npts';Npts(end:-1:1)'],[err_max;err_min(end:-1:1)],h1.Color);
end
h.FaceAlpha = 0.15;

% for i = ind
%     loglog(Npts,err_reg(:,i).^2,'--','linewidth',1,'Color',0.75*ones(1,3));
%     loglog(Npts,err_emb(:,i).^2,'-','linewidth',1,'Color',0.75*ones(1,3));
% end

% reference
ylabel '$$d^2_{Gr(p,n)}\left(s,\hat{s}\right)$$'
xlabel 'number of landmarks $$(n)$$';
legend('$$O(n^{-4})$$','interpreter','latex','Location','southwest')

% compare curvatures
% figure;
% % embedding representation
% h1 = semilogy(emb.L/max(emb.L),emb.curv,'linewidth',2); hold on; grid on;
% % Cartesian spline
% h2 = semilogy(tru_t,curv2,'--','linewidth',2);
% xlabel '$$t$$'; ylabel '$$ \hat{\kappa}(t)$$';

% compare shapes
% figure;
% plot(reg_pts(:,1),reg_pts(:,2),'--','linewidth',2); hold on; axis equal;
% plot(emb_pts(:,1),emb_pts(:,2),'--','linewidth',2);
% l = linspace(0,1,5000); [coordU, coordL] = cst_airfoil(l',dv(1:m/2)',dv(m/2+1:end)',0);
% bnd = [coordL'; coordU(:,end:-1:1)']; x0 = [mean(bnd(:,1)),mean(bnd(:,2))]; bnd = bnd - x0;
% plot(bnd(:,1),bnd(:,2),'k','linewidth',2); 
% plot(P0(:,1),P0(:,2),'ko')

% plot embedding representation approximation
fig = figure; 
bndplot(fig,emb.pts(:,1),emb.pts(:,2),log(emb.curv)); hold on;
plot(emb_pts(:,1),emb_pts(:,2),'k','linewidth',2); colorbar off;
% fig.CurrentAxes.Visible = 'off';

