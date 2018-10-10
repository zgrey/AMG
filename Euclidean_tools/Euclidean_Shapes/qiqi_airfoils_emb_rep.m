% post-process QiQi airfoil database to check radial convexity of shapes
clc; clearvars; close all;
addpath ../;
%% Things to modify
% plot flag to visualize and print the embedding
plt_flg = 0;
% number of resampled landmarks
N = 500;

%% read airfoil
dirstruct  = dir('~/airfoil_coordinate_database/coordinates/');
Naf = size(dirstruct,1)-2;
ind_emb = ones(Naf,1); nonconvcount = ones(Naf,1);
% precondition alpha and theta matrices for KLE
aa = zeros(Naf,N); th = zeros(Naf,N);
for i=1:Naf
clc; fprintf('%0.2f%% complete...\n',i/(Naf)*100);
%% read nominal points from file
P0 = importairfoil([dirstruct(2+i).folder,'/',dirstruct(2+i).name]);
% number of points
n = size(P0,1);

%% affine transformations (i.e., p  = M*p0   + a)
% shift and scale to [-1,1] box
% pu = max(P0,[],1)'; pl = min(P0,[],1)';
% M = diag(2./(pu-pl)); 
% a = -(M*pl + ones(2,1)); 

% shift and scale using landmark-affine standardization (Bryner, 2D affine and projective spaces)
C_x = mean(P0,1)'; Sig_x = (P0 - repmat(C_x',n,1))'*(P0 - repmat(C_x',n,1));
M = chol(Sig_x) \ eye(2); a = -M*C_x;
% affine inverse (i.e., x0 = Minv*x + ainv in [xl,xu]_m)
Minv = chol(Sig_x); ainv = C_x; 

% shift and scale to center of mass
% pu = max(P0,[],1)'; pl = min(P0,[],1)';
% M = diag(2./(pu-pl));
% a = -M*mean(P0,1)';

% transform points
P = (P0*M + repmat(a',n,1));

% repeat first point to close shape
P = [P;P(1,:)];

%% embedding
% compute discrete lengths
t = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 + ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
% compute unique angles using atan2
% s0 = atan2(P(:,2),P(:,1)); s0 = unwrap(s0);
% resolve discontinuity from atan2 by summing angular discrepancies
ds0 = [atan2(P(1,2),P(1,1));
       atan2(P(1:end-1,1).*P(2:end,2) - P(1:end-1,2).*P(2:end,1),...
             P(1:end-1,1).*P(2:end,1) + P(1:end-1,2).*P(2:end,2))];
s0 = cumsum(ds0);

% scale
s0 = s0/(2*pi) + 1/2;

% circular
a = 1; b = 1;
nu_emb = [b*cos(2*pi*(s0-1/2)) , a*sin(2*pi*(s0-1/2))];
% compute inner product
a0 = P(:,1).*nu_emb(:,1) + P(:,2).*nu_emb(:,2);
% take only unique points
[t,ia] = unique(t,'stable');
a0 = a0(ia); s0 = s0(ia);
% construct radial embedding
[ind,S,nu,S_emb,a_spl,s_spl,s,alpha,t]= radial_emb3_sort(s0,a0,t,N);
% count non-radially convex airfoils
ind_emb(i) = sum(ind ~= 0);
% save uniform discretization of inner product and angular distributions
aa(i,:) = ppval(a_spl,linspace(0,1,N));
th(i,:) = ppval(s_spl,linspace(0,1,N));

%% convex check
k = convhull(P(:,1),P(:,2));
nonconvcount(i) = n - length(k);

%% transform back to original coordinates
% shape landmarks
S_emb0 = S_emb*Minv + repmat(ainv',N,1);
% plot continuous approximation of normals
nu0 = nu*[0 -1; 1 0]*Minv*[0 1; -1 0];
% compute unit normals
nu0 = nu0./sqrt(nu0(:,1).^2 + nu0(:,2).^2);

%% visualize
if plt_flg == 1
figure(1);
% plot embedding representation
subplot(1,2,1), h2 = plot(S_emb(:,1),S_emb(:,2),'linewidth',2); axis equal; hold on;
% plot approximate unit normals
subplot(1,2,1), h6 = quiver(S_emb(:,1),S_emb(:,2),nu(:,1),nu(:,2),'color',h2(1).Color);
% plot nominal points
subplot(1,2,1), h1 = scatter(P(:,1),P(:,2),20,'k');
% plot embedding at nominal points (check interpolation)
subplot(1,2,1), h9 = plot(S(:,1),S(:,2),'.','MarkerEdgeColor',h2(1).Color,'MarkerSize',10);
% plot non-radially convex points
subplot(1,2,1), h7 = plot(P(ind,1),P(ind,2),'r.','MarkerSize',10);
% plot convex hull
subplot(1,2,1), h8 = plot(P(k,1),P(k,2),'--','color',0.75*ones(1,3));
% plot 3-dim. embedding
subplot(1,2,2), h3 = plot3(s,alpha,t,'ko',ppval(s_spl,linspace(0,1,N)),ppval(a_spl,linspace(0,1,N)),linspace(0,1,N)); grid on; hold on;
% plot non-radially convex points in embedding (non-monotone angles)
subplot(1,2,2), h4 = scatter3(s(ind(ind ~= 0)),alpha(ind(ind ~= 0)),t(ind(ind ~= 0)),'filled','r');
xlabel 's(t)'; ylabel '\alpha(t)'; zlabel 't';
h5 = annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Directory index i=',num2str(i),': ',dirstruct(2+i).name], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
% print figure
print(['./figs_emb3/',num2str(i),'_',dirstruct(2+i).name,'.png'],'-dpng')
% clear plot
delete([h1;h2;h3;h4;h5;h6;h7;h8;h9]); subplot(1,2,1), reset(gca); subplot(1,2,2), reset(gca);

fig = figure(2); hold on; axis equal;
% plot continuous approximation of shape
h1 = plot(S_emb0(:,1),S_emb0(:,2),'linewidth',2); axis equal; hold on;
h2 = quiver(S_emb0(:,1),S_emb0(:,2),nu0(:,1),nu0(:,2),'color',h1.Color);
% plot original airfoil
PP0 = P*Minv + repmat(ainv',n+1,1);
h3 = scatter(PP0(:,1),PP0(:,2),15,'k');
fig.CurrentAxes.Visible = 'off';

h4 = annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Directory index i=',num2str(i),': ',dirstruct(2+i).name], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
% print figure
print(['./figs_airfoils/',num2str(i),'_',dirstruct(2+i).name,'.png'],'-dpng')
% clear plot
delete([h1;h2;h3;h4]); reset(gca);
end
end
fprintf('%f%% (%i / %i) are radially convex\n',(Naf - length(find(ind_emb)))/Naf*100,(Naf - length(find(ind_emb))),Naf)
fprintf('%f%% (%i / %i) are convex\n',(Naf - length(find(nonconvcount)))/Naf*100,(Naf - length(find(nonconvcount))),Naf)

%% construct Karhunen-Loeve parameterization
nrnd = 100;
% Karhunen-Loeve Expansion of inner products
[Ua,Da,Va,Arnd,a_mu,a_cdf,a_z] = KLE(aa,nrnd);
% Karhunen-Loeve Expansion of angle
[Uth,Dth,Vth,THrnd,th_mu,th_cdf,th_z] = KLE(th,nrnd);

%% Visualize KLE
% random evaluations
figure;
% plot 3-dim. embedding
plot3(th_mu,a_mu,linspace(0,1,N)','k',ppval(s_spl,linspace(0,1,N)'),ppval(a_spl,linspace(0,1,N)'),linspace(0,1,N)'); grid on; hold on;
indx = randi(nrnd,10,1); plot3(THrnd(:,indx),Arnd(:,indx),linspace(0,1,N)','color',0.25*ones(1,3));

% empirical distribution of inner product KL coefficients
figure; hold on; grid on;
for i=1:N, plot(a_z(:,i),a_cdf(:,i),'color',0.65*ones(1,3)); end
plot(linspace(-5,5,1000),normcdf(linspace(-5,5,1000)),'k','linewidth',2); hold on;
xlabel '\alpha'; title 'Standardized Empirical CDFs of Karhunen-Loeve Coefficients';

% empirical distribution of KL coefficients
figure; hold on; grid on;
for i=1:N, plot(th_z(:,i),th_cdf(:,i),'color',0.65*ones(1,3)); end
plot(linspace(-5,5,1000),normcdf(linspace(-5,5,1000)),'k','linewidth',2); hold on;
xlabel '\theta';  title 'Standardized Empirical CDFs of Karhunen-Loeve Coefficients';

% evaluate random shape
a_mu = Arnd(:,50); th_mu = THrnd(:,50);
S_mu =  [a_mu.*cos(2*pi*(th_mu-1/2)),...
         a_mu.*sin(2*pi*(th_mu-1/2))];
% compute continuous normal vector approximation
nu_mu  = a_mu.*([sin(2*pi*(th_mu-1/2)), -cos(2*pi*(th_mu-1/2))]) +...
         2*pi*th_mu.*[a_mu.*cos(2*pi*(th_mu-1/2)), a_mu.*sin(2*pi*(th_mu-1/2))];        
% transform back to original coordinates
S_mu0 = S_mu*Minv + repmat(ainv',N,1);
% plot continuous approximation of normals
nu0_mu = nu_mu*[0 -1; 1 0]*Minv*[0 1; -1 0];
% compute unit normals
nu0_mu = nu0_mu./sqrt(nu0_mu(:,1).^2 + nu0_mu(:,2).^2);
% plot the resulting average shape
figure; plot(S_mu0(:,1),S_mu0(:,2),'k','LineWidth',2); axis equal;