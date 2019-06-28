% convergence study
clc; close all; clearvars; rng(42);
addpath ../../../ASAP/AUTO

% number of resampled points (also used as total number of points in Riemann sums of error functionals)
NN = 10^4;
% maximum order of magnitude increase for convergence study
O = 5;

%% test shapes
% circle test case
% Pshp =@(t) [cos(t*2*pi),sin(t*2*pi)];
% shp_curv =@(t) ones(length(t),1);

% ellipse test case
% G = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
% scl = 2; Pshp =@(t) [scl*cos(2*pi*t),sin(2*pi*t)]*G'; 
% shp_curv =@(t) scl./( scl^2*sin(2*pi*t).^2 + cos(2*pi*t).^2).^(3/2);

% peanut test case
% G = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
% scl = 0.5; Pshp =@(t) [scl*(cos(4*pi*t)+2).*cos(2*pi*t), (cos(4*pi*t)+2).*sin(2*pi*t)]*G';
% shp_curv = @(t) 24*sqrt(2)*scl*abs(-7-8*cos(4*pi*t)+cos(8*pi*t))./...
%                ( 18 + 34*scl^2 + (27+5*scl^2)*cos(4*pi*t) - 6*(5*scl^2-3)*cos(8*pi*t) + 9*cos(12*pi*t) - 9*scl^2*cos(12*pi*t) ).^(3/2);

% parameterized star test case
Nshp = 99;
err_emb = ones(O,Nshp); err_reg = ones(O,Nshp);
for ii = 0:Nshp
% G = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
G = eye(2);
scl = 1; w = ii; Pshp =@(t) [scl*(cos(w*2*pi*t)+2).*cos(2*pi*t), (cos(w*2*pi*t)+2).*sin(2*pi*t)]*G';
shp_curv = @(t) abs(scl*(3*(12+w^2) + 4*(8+w^2)*cos(pi*t*w) - (w^2-4)*cos(2*pi*t*w)))./...
                (sin(2*pi*t).^2.*(16*scl^2 + 16*scl^2*cos(pi*t*w) + 4*scl^2*cos(pi*t*w).^2 + w^2*sin(pi*t*w).^2)...
                + cos(2*pi*t).^2.*(16 + 16*cos(pi*t*w) + 4*cos(pi*t*w).^2 + scl^2*w^2*sin(pi*t*w).^2)...
                + (scl^2-1)*w*sin(4*pi*t).*(4*sin(pi*t*w) + sin(2*pi*t*w))...
                ).^(3/2);

% preallocate 
err1 = ones(O,3); err2 = ones(O,3); Npts = ones(O,1);
for i=1:O
    % generate uniform spacing of landmarks from shape
    % repeats last point
%     t0 = 0:10^(-i)*2.5:1; P0 = Pshp(t0'); Npts(i) = length(t0);
    t0 = 0:10^(-i):1; P0 = Pshp(t0'); Npts(i) = length(t0);
    
    %% 3D circular embedding
    % get default options
    options = embOptions([]);
    % compute embedding representation
    emb = embrep3(P0,NN,'uni',options);
    
    %% regular spline comparison
%     P = unique_points(P0); P = [P; P(1,:)];
    % if data has the known repeated last point
    P = P0; rmv = 0;
    % if data does not have a repeated point (repeat the first for periodicity)
%     P = [P0; P0(1,:)]; rmv = 1;
    
    % discrete length scales for splines
    l0 = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 + ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
    ub = max(l0); lb = 0;
    
    % scale to [0,1]
    t = l0/ub;
    x_spl = csape(t,P(:,1),'periodic');
    y_spl = csape(t,P(:,2),'periodic');
    
    % recompute points at refined landmarks
    tt = linspace(0,1,NN)';
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
 
    % reevaluate at points for error computation
    % convert between length scale measures
    t_spl = pchip(emb.nom.t,t0);
    treval = ppval(t_spl,emb.t);
    emb_curv = emb.curv; emb_pts = emb.pts;
    reg_curv = curv2; reg_pts = P;
    % recompute curvature and points
    tru_crv = shp_curv(treval);
    PP = Pshp(treval);
    % compute error
    P0TF = affine_trans(PP(1:end-1,:),'LA');
    embTF = affine_trans(emb.pts(1:end-1,:),'LA');
    regTF = affine_trans(reg_pts(1:end-1,:),'LA');
    % NOTE NN-2 IS THE NUMBER OF UNIQUE POINTS (SCALING DIFFERENTLY CHANGES THE ERROR METRIC)
    err1(i,:) = [1/NN*sqrt(sum( (tru_crv - emb_curv).^2, 1)), 1/NN*sum(sqrt(sum((PP - emb_pts).^2,2))),...
        dGr_np(1/sqrt(NN-2)*P0TF,1/sqrt(NN-2)*embTF)];
    err2(i,:) = [1/NN*sqrt(sum( (tru_crv - reg_curv).^2, 1)), 1/NN*sum(sqrt(sum((PP - reg_pts).^2,2))),...
        dGr_np(1/sqrt(NN-2)*P0TF,1/sqrt(NN-2)*regTF)];
    
    clc; disp(['AGGREGATING DATA...',num2str(i/O*100),'% Complete']);
    disp([num2str(ii+1),'/',num2str(Nshp+1),' Shapes Complete']);
end
    err_emb(:,ii+1) = err1(:,3);
    err_reg(:,ii+1) = err2(:,3);
end
%% convergence plots
set(0,'defaulttextInterpreter','latex','DefaultAxesFontSize',24)
figure; 
% loglog(Npts,Npts.^(-2),'k','linewidth',2); grid on; hold on; 
% curvature errors
% h2 = loglog(Npts,err_reg(:,1),'-d','linewidth',2);
% h1 = loglog(Npts,err_emb(:,1),'-o','linewidth',2); 
% % landmark errors
% loglog(Npts,err_emb(:,2),'--o','linewidth',2,'color',h1.Color); 
% loglog(Npts,err_reg(:,2),'--d','linewidth',2,'color',h2.Color);
% Grassmannian distance
% loglog(Npts,err2(:,3).^2,'-d','linewidth',2,'color',h2.Color);
% loglog(Npts,err1(:,3).^2,'-o','linewidth',2,'color',h1.Color); 
% ylabel '$$\Vert \hat{\kappa}_{c}(t_i) - \kappa_{c}(t_i) \Vert_2$$';

h1 = loglog(Npts,mean(err_reg,2).^2,'--d','linewidth',2,'markersize',8); grid on; hold on; 
h2 = loglog(Npts,mean(err_emb,2).^2,'-o','linewidth',2,'markersize',8); 
% plot some random shape convergence
% err_sub = err_emb.^2; ind = err_sub < eps; err_sub(ind) = 1/2*eps;
% loglog(Npts,err_sub(:,1:10:Nshp+1),'-o','linewidth',1,'Color',0.75*ones(1,3));
% err_sub = err_reg.^2; ind = err_sub < eps; err_sub(ind) = 1/2*eps;
% loglog(Npts,err_sub(:,1:10:Nshp+1),'--d','linewidth',1,'Color',0.75*ones(1,3));

err_max = max(err_reg,[],2).^2; err_min = min(err_reg,[],2).^2;
% replace zeros with significantly reduced error for plots
ind1 = err_max == 0; err_max(ind1) = 1/2*eps;
ind2 = err_min == 0; err_min(ind2) = 1/2*eps;
h = fill([Npts(1:end);Npts(end:-1:1)],[err_max(1:end);err_min(end:-1:1)],h1.Color);
scatter([Npts(ind1);Npts(ind2)],[err_max(ind1),err_min(ind2)],'ko','linewidth',2)
h.FaceAlpha = 0.15;
err_max = max(err_emb,[],2).^2; err_min = min(err_emb,[],2).^2;
% replace zeros with significantly reduced error for plots
ind1 = err_max == 0; err_max(ind1) = 1/2*eps;
ind2 = err_min == 0; err_min(ind2) = 1/2*eps;
h = fill([Npts(1:end);Npts(end:-1:1)],[err_max(1:end);err_min(end:-1:1)],h2.Color);
scatter([Npts(ind1);Npts(ind2)],[err_max(ind1),err_min(ind2)],'ko','linewidth',2)
h.FaceAlpha = 0.15;

ylabel '$$d^2_{Gr(p,n)}\left(s,\hat{s}\right)$$'
xlabel 'number of landmarks $$(n)$$';
% legend('$$O(n^{-2})$$','interpreter','latex','Location','southwest')
chH = get(gca,'Children');
set(gca,'Children',[chH(end-1);chH(end);chH(1:end-2)])

% compare curvatures
% figure;
% % known curvature
% if exist('tru_crv','var')
% semilogy(linspace(0,1,NN)',shp_curv(linspace(0,1,NN)'),'k','linewidth',2); hold on;
% semilogy(treval,tru_crv,'ko');
% end
% % embedding representation
% h1 = semilogy(treval,emb_curv,'*-','linewidth',2); hold on; grid on;
% % semilogy(emb.nom.th/(2*pi),emb.nom.curv,'*','color',h1.Color);
% % Cartesian spline
% % h2 = semilogy(treval,reg_curv,'d-','linewidth',2);
% % semilogy(emb.nom.th/(2*pi),nom_curv2,'o','color',h2.Color)
% xlabel '$$t$$'; ylabel '$$ \hat{\kappa}(t)$$';

% compare shapes
% figure; skp = floor(0.01*NN);
% plot(PP(1:skp:end,1),PP(1:skp:end,2),'k*','markersize',10); hold on; axis equal;
% plot(emb_pts(1:skp:end,1),emb_pts(1:skp:end,2),'*--','linewidth',2)
% plot(reg_pts(1:skp:end,1),reg_pts(1:skp:end,2),'d--','linewidth',2,'markersize',10)

% plot embedding representation approximation
fig = figure; 
bndplot(fig,emb_pts(:,1),emb_pts(:,2),log(emb_curv)); hold on;
plot(emb_pts(:,1),emb_pts(:,2),'k','linewidth',2)
fig.CurrentAxes.Visible = 'off';

