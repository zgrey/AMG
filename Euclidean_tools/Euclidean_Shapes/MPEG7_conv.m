% convergence study
clc; close all; clearvars; rng(42);

%% image processing test shapes

dirstruct = dir('./MPEG7dataset/original/');

Nshp = size(dirstruct,1)-2;
for ii = 1:Nshp
% Rim = imread('./ralphie.png'); Nshp = 1;
% Rim = imread('./MPEG7dataset/original/octopus-16.gif'); Nshp = 1;
% Rim = imread('./MPEG7dataset/original/bat-1.gif'); Nshp = 1;
% Rim = imread('./MPEG7dataset/original/dog-1.gif'); Nshp = 1;
% Rim = imread('./MPEG7dataset/original/beetle-6.gif'); Nshp = 1;

% read image
Rim = imread([dirstruct(2+ii).folder,'/',dirstruct(2+ii).name]);
BW = imbinarize(Rim(:,:,1)); Rbnd = bwboundaries(BW);
% look for image boundary (assumed to have the most landmarks)
if size(Rbnd,1) > 1
    k = 1;
    for iii = 2:size(Rbnd,1)
        if size(Rbnd{iii},1) > size(Rbnd{k},1)
            k = iii;
        end
    end
else
    k = 1;
end
bnd = Rbnd{k}*[0 -1;1 0]; % for MPEG7 database
% flip to be counter-clockwise
bnd = bnd(end:-1:1,:);

[bnd,Ntru] = unique_points(bnd); bnd = [bnd; bnd(1,:)];
tru_t = cumsum([0; sqrt(( bnd(2:end,1) - bnd(1:end-1,1) ).^2 + ( bnd(2:end,2) - bnd(1:end-1,2) ).^2)],1);
tru_t = tru_t/max(tru_t); 

%% shift and scale
x0 = [mean(bnd(:,1)),mean(bnd(:,2))]; bnd = bnd - x0;

Nmax = size(bnd,1);
skp = 2.^(1:floor(log(Nmax)/log(2))-1); skp = skp(end:-1:1);

%% Convergence study
for i=1:length(skp)
    % select subset (ensure point set is from boundary of a simply-connected domain)
    P0 = bnd(1:skp(i):end,:); 
    P0 = unique_points(P0); P0 = [P0; P0(1,:)]; Npts(i,ii) = size(P0,1);
    P = P0;
    
    %% 3D circular embedding
    % get default options
    options = embOptions([]); 
%     options.AffineTrans = 'none';
    % compute embedding representation
    emb = embrep3(P0,Nmax,'nom',options);
    
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

%% Fit average errors
NN = reshape(Npts,size(Npts,1)*Nshp,1);
zero_ind = NN == 0;
NN = NN(~zero_ind);
err1 = reshape(err_emb,size(err_emb,1)*Nshp,1);
ferr_emb = err1(~zero_ind).^2;
err2 = reshape(err_reg,size(err_reg,1)*Nshp,1);
ferr_reg = err2(~zero_ind).^2;

coefs1 = [ones(size(NN,1),1), log(NN)] \ log(ferr_emb);
coefs2 = [ones(size(NN,1),1), log(NN)] \ log(ferr_reg);
res1 = zeros(size(err_emb,1),size(err_emb,2));
for ii=1:Nshp
res1(:,ii) = log(err_emb(:,ii).^2) - ([ones(size(Npts(:,ii),1),1), log(Npts(:,ii))]*coefs1);
end
ind = res1 == -Inf;
res1(ind) = 0;
%% convergence plots
set(0,'defaulttextInterpreter','latex','DefaultAxesFontSize',24)
figure; NNN1 = exp(linspace(log(min(NN)),log(max(NN)),4));
NNN2 = exp(linspace(log(min(NN)),log(max(NN)),3));
loglog(NNN1,NNN1.^(-2),'k','linewidth',2); hold on; grid on;


% Average Cartesian distances (errors)
h1 = loglog(NNN2,exp(coefs2(1) + coefs2(2)*log(NNN2)),'--d','linewidth',2,'markersize',8);
% Average embedding distances (errors)
h2 = loglog(NNN1, exp(coefs1(1) + coefs1(2)*log(NNN1)),'-o','linewidth',2,'markersize',8);

% plot some shape convergence
loglog(Npts(:,1:100:Nshp),err_emb(:,1:100:Nshp).^2,'-','linewidth',0.5,'Color',0.75*ones(1,3));
loglog(Npts(:,1:100:Nshp),err_reg(:,1:100:Nshp).^2,'--d','linewidth',0.5,'Color',0.75*ones(1,3));

[~,ind2] = min(sum(res1,1)); [~,ind1] = max(sum(res1,1));
err_max = err_emb(:,ind1).^2; ind_max = err_max ~= 0; err_max = err_max(ind_max); Nub = Npts(ind_max,ind1);
err_min = err_emb(:,ind2).^2; ind_min = err_min ~= 0; err_min = err_min(ind_min); Nlb = Npts(ind_min,ind2);
h = fill([Nub; Nlb(end:-1:1)],[err_max; err_min(end:-1:1)],h2.Color);
h.FaceAlpha = 0.15;

err_max = err_reg(:,ind1).^2; ind_max = err_max ~= 0; err_max = err_max(ind_max); Nub = Npts(ind_max,ind1);
err_min = err_reg(:,ind2).^2; ind_min = err_min ~= 0; err_min = err_min(ind_min); Nlb = Npts(ind_min,ind2);
h = fill([Nub; Nlb(end:-1:1)],[err_max; err_min(end:-1:1)],h1.Color);
h.FaceAlpha = 0.15;

ylabel '$$d^2_{Gr(p,n)}\left(s,\hat{s}\right)$$'
xlabel 'number of landmarks $$(n)$$';
legend('$$O(n^{-2})$$','interpreter','latex','Location','southwest')
chH = get(gca,'Children');
set(gca,'Children',[chH(end-1);chH(end);chH(1:end-2)])

% compare curvatures
% figure;
% % embedding representation
% h1 = semilogy(emb.L/max(emb.L),emb.curv,'linewidth',2); hold on; grid on;
% % Cartesian spline
% h2 = semilogy(tru_t,curv2,'--','linewidth',2);
% xlabel '$$t$$'; ylabel '$$ \hat{\kappa}(t)$$';

% compare shapes
% figure;
% plot(bnd(:,1),bnd(:,2),'linewidth',2); hold on; axis equal;
% plot(bnd(:,1),bnd(:,2),'w.');
% plot(P0(:,1),P0(:,2),'ko')
% plot(emb_pts(:,1),emb_pts(:,2),'-.','linewidth',2)
% plot(reg_pts(:,1),reg_pts(:,2),'--','linewidth',2)

% plot embedding representation approximation
% fig = figure; 
% bndplot(fig,emb_pts(:,1),emb_pts(:,2),log(emb.curv)); hold on; axis equal;
% plot(emb.pts(:,1),emb.pts(:,2),'w.')

