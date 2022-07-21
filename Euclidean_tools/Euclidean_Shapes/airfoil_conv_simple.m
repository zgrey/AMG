% convergence study
clc; close all; clearvars; rng(41);

%% airfoil study
addpath ../../../ASAP/AUTO/;
% optimize to find nominal values
m = 18; nl = 5001; 
% uniform over parameterization domain (not uniform over shape!)
% l = linspace(0,1,nl);
% an "intuitive" shape sampling (clustered at LE and TE)
l = (cos(linspace(0,pi,nl)) + 1)/2;
% an "unintuitive" shape sampling (clustered at mid-chord)
% l = (sin(linspace(pi/2,-pi/2,nl)) + 1)/2;

dv0 = 0.45*[-ones(m/2,1); ones(m/2,1)];
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
    % sample random shape
    dv = bsxfun(@plus,lb0,bsxfun(@times,ub0-lb0,0.5*(X(ii,:)+1)));
    [coordU, coordL] = cst_airfoil(l',dv(1:m/2)',dv(m/2+1:end)',0);
    bnd = [coordL'; coordU(:,end:-1:1)'];
    [bnd,Ntru] = unique_points(bnd); bnd = [bnd; bnd(1,:)];
    tru_t = cumsum([0; sqrt(( bnd(2:end,1) - bnd(1:end-1,1) ).^2 +...
                    ( bnd(2:end,2) - bnd(1:end-1,2) ).^2)],1);
    tru_t = tru_t/max(tru_t); 

    %% shift and scale
    % x0 = [mean(bnd(:,1)),mean(bnd(:,2))]; bnd = bnd - x0;

    Nmax = size(bnd,1);
    skp = 2.^(0:floor(log(Nmax)/log(2))-2); skp = skp(end:-1:1);

    %% Convergence study
    for i=1:length(skp)-1
        % generate new (sparse) shape
        N_sparse = length(1:skp(i):(Nmax-1)/2 +1);
        l_sparse = (cos(linspace(0,pi,N_sparse)) + 1)/2;
        [coordU_sparse, coordL_sparse] = cst_airfoil(l_sparse',dv(1:m/2)',dv(m/2+1:end)',0);
        bnd_sparse = [coordL_sparse'; coordU_sparse(:,end:-1:1)'];
        
        % assign unique sparse landmarks for analysis (double check count)
        P0 = unique_points(bnd_sparse); P0 = [P0; P0(1,:)]; Npts(i) = size(P0,1)-1;
        P = P0;

        %% 3D circular embedding
        % get default options
        options = embOptions([]); 
        % compute embedding representation
        emb = embrep3(P0,250,'nom',options);

        %% regular spline comparison
        % discrete length scales for splines
        t0 = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 +...
                     ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
        hbar(i) = max(sqrt(sum(diff(P,1,1).^2,2)));
        ub = max(t0); lb = 0;

        % scale to [0,1]
        t = t0/ub;
        x_spl = csape(t,P(:,1),'periodic');
        y_spl = csape(t,P(:,2),'periodic');

        % recompute points at refined "true" landmarks
        P = [ppval(x_spl,tru_t), ppval(y_spl,tru_t)];

        %% reevaluate at "consistent" points for error computation
        % convert between length scale measures
        treval = pchip(t,emb.TF.t,tru_t);
        
        emb_pts = (ppval(emb.alph.spl,treval).*[cos(ppval(emb.th.spl,treval)),...
                   sin(ppval(emb.th.spl,treval))] ...
                   - repmat(emb.TF.b',length(tru_t),1))*emb.TF.Minv';
        emb_pts_misreg = (ppval(emb.alph.spl,t).*[cos(ppval(emb.th.spl,t)),...
                   sin(ppval(emb.th.spl,t))] ...
                   - repmat(emb.TF.b',length(t),1))*emb.TF.Minv';
        emb_pts_reg = (ppval(emb.alph.spl,emb.TF.t).*[cos(ppval(emb.th.spl,emb.TF.t)),...
                   sin(ppval(emb.th.spl,emb.TF.t))] ...
                   - repmat(emb.TF.b',length(emb.TF.t),1))*emb.TF.Minv';
        reg_pts = P;

        %% compute shape disctance over the Grassmannian
        % transform the "true" refined shape
        bndTF = affine_trans(bnd(1:end-1,:),'LA');
        % transform the refined embedding representation
        embTF = affine_trans(emb_pts(1:end-1,:),'LA');
        dGr_emb(i,ii) = dGr_np(1/sqrt(Ntru-1)*bndTF,1/sqrt(Ntru-1)*embTF);
        err_emb(i,ii) = max(sqrt(sum((bnd - emb_pts).^2,2)));
%         err_emb(i,ii) = norm(bnd - emb_pts,"inf");
%         err_emb(i,ii) = max(sqrt(sum((bndTF - embTF).^2,2)));
%         err_emb(i,ii) = dGr_emb(i,ii);
        % transform the refined regular (planar spline) representation
        regTF = affine_trans(reg_pts(1:end-1,:),'LA');
        dGr_reg(i,ii) = dGr_np(1/sqrt(Ntru-1)*bndTF,1/sqrt(Ntru-1)*regTF);
        err_reg(i,ii) = max(sqrt(sum((bnd - reg_pts).^2,2)));
%         err_reg(i,ii) = norm(bnd - reg_pts,"inf");
%         err_reg(i,ii) = max(sqrt(sum((bndTF - regTF).^2,2)));
%         err_reg(i,ii) = dGr_reg(i,ii);

        clc; disp(['AGGREGATING DATA...',num2str(i/length(skp)*100),'% Complete']);
        disp([num2str(ii),'/',num2str(Nshp),' Shapes Complete']);
    end
end

%% convergence plots
set(0,'defaulttextInterpreter','latex','DefaultAxesFontSize',24)
figure; 

% loglog(Npts,Npts.^(-4),'k--','linewidth',2); hold on;
loglog(hbar,hbar.^(4),'k--','linewidth',2); hold on;
loglog(hbar,hbar.^(2),'k--','linewidth',2); hold on;

% Cartesian distances (errors)
% h1 = loglog(Npts,mean(dGr_reg.^2,2),'--d','linewidth',2,'markersize',8);  grid on; hold on;
% h3 = loglog(Npts,mean(err_reg,2),'--d','linewidth',2,'markersize',8);  grid on; hold on;
h1 = loglog(hbar,mean(err_reg,2),'--d','linewidth',2,'markersize',8);  grid on; hold on;
% embedding distances (errors)
% h2 = loglog(Npts,mean(err_emb,2),'-o','linewidth',2,'markersize',8);
h2 = loglog(hbar,mean(err_emb,2),'-o','linewidth',2,'markersize',8);

err_max = max(err_reg,[],2); err_min = min(err_reg,[],2);
if min(err_min) == 0
%     h = fill([Npts';Npts(end:-1:1)'],[err_max(1:end);1/2*eps;err_min(end-1:-1:1)],h1.Color);
%     scatter(Npts(end),1/2*eps,'ko','linewidth',2)
%       h = fill([hbar';hbar(end:-1:1)'],[err_max(1:end);1/2*eps;err_min(end-1:-1:1)],h1.Color);
%       scatter(hbar(end),1/2*eps,'ko','linewidth',2)
else
%     h = fill([Npts';Npts(end:-1:1)'],[err_max;err_min(end:-1:1)],h1.Color);
%     h = fill([hbar';Npts(end:-1:1)'],[err_max;err_min(end:-1:1)],h1.Color);
end
% h.FaceAlpha = 0.15;
err_max = max(err_emb,[],2); err_min = min(err_emb,[],2);
if min(err_min) == 0
%     h = fill([Npts';Npts(end:-1:1)'],[err_max(1:end);1/2*eps;err_min(end-1:-1:1)],h2.Color);
%     scatter(Npts(end),1/2*eps,'ko','linewidth',2)
%       h = fill([hbar';hbar(end:-1:1)'],[err_max(1:end);1/2*eps;err_min(end-1:-1:1)],h2.Color);
%       scatter(hbar(end),1/2*eps,'ko','linewidth',2)
else
%     h = fill([Npts';Npts(end:-1:1)'],[err_max;err_min(end:-1:1)],h1.Color);
%       h = fill([hbar';Npts(end:-1:1)'],[err_max;err_min(end:-1:1)],h1.Color);
end
% h.FaceAlpha = 0.15;

% ylabel '$$d^2_{\mathcal{G}(1e4,2)}\left([\tilde{X}],[\tilde{X}^*]\right)$$'
ylabel '$$\Vert \tilde{X} - \tilde{X}^*\Vert_{\infty,2}$$'
% xlabel 'number of landmarks $$(n_c)$$';
xlabel 'gauge $$(\bar{h})$$';
% legend('$$O(n_c^{-4})$$','interpreter','latex','Location','southwest')
legend('$$O(\bar{h}^{4})$$','interpreter','latex','Location','southwest')
