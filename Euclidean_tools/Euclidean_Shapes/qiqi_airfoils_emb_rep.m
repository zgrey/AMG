% post-process QiQi airfoil database to check radial convexity of shapes
clc; clearvars; close all;
% read airfoil
dirstruct  = dir('~/airfoil_coordinate_database/coordinates/');
Naf = size(dirstruct,1)-2;
ind_emb = ones(Naf,1); nonconvcount = ones(Naf,1);
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
M = chol(Sig_x) \ eye(2); 
a = -M*C_x;

% shift and scale to center of mass
% pu = max(P0,[],1)'; pl = min(P0,[],1)';
% M = diag(2./(pu-pl));
% a = -M*mean(P0,1)';

% affine inverse (i.e., x0 = Minv*x + ainv in [xl,xu]_m)
Minv = chol(Sig_x); ainv = C_x;  

% transform points
P = (P0*M + repmat(a',n,1));

% repeat first point to close shape
P = [P;P(1,:)];

%% embedding
% compute discrete lengths
t = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 + ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
% compute unique angles using atan2
% resolve discontinuity from atan2 by summing angular discrepancies
% s0 = atan2(P(:,2),P(:,1)); s0 = unwrap(s0);
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
% sort points based on embedding
N = 500;
% [ind,S,nu,S_rc,pp,s,alp]= radial_emb_sort(s0,a0,N,10,0,a,b);
[ind,S,nu,S_emb,pp,ss,s,alph,t]= radial_emb3_sort(s0,a0,t,N);
ind_emb(i) = sum(ind ~= 0);

%% convex check
k = convhull(P(:,1),P(:,2));
nonconvcount(i) = n - length(k);

%% visualize
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
subplot(1,2,2), h3 = plot3(s,alph,t,'ko',ppval(ss,linspace(0,1,N)),ppval(pp,linspace(0,1,N)),linspace(0,1,N)); grid on; hold on;
% plot non-radially convex points in embedding (non-monotone angles)
subplot(1,2,2), h4 = scatter3(s(ind(ind ~= 0)),alph(ind(ind ~= 0)),t(ind(ind ~= 0)),'filled','r');
xlabel 's(t)'; ylabel '\alpha(t)'; zlabel 't';
h5 = annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Directory index i=',num2str(i),': ',dirstruct(2+i).name], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
% print figure
print(['./figs_emb3/',num2str(i),'_',dirstruct(2+i).name,'.png'],'-dpng')
% clear plot
delete([h1;h2;h3;h4;h5;h6;h7;h8;h9]); subplot(1,2,1), reset(gca); subplot(1,2,2), reset(gca);

%% transform back to original coordinates
% shape landmarks
S_emb0 = S_emb*Minv + repmat(ainv',N,1);
% plot continuous approximation of normals
nu0 = nu*[0 -1; 1 0]*Minv*[0 1; -1 0];
% compute unit normals
nu0 = nu0./sqrt(nu0(:,1).^2 + nu0(:,2).^2);

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
fprintf('%f%% (%i / %i) are radially convex\n',(Naf - length(find(ind_emb)))/Naf*100,(Naf - length(find(ind_emb))),Naf)
fprintf('%f%% (%i / %i) are convex\n',(Naf - length(find(nonconvcount)))/Naf*100,(Naf - length(find(nonconvcount))),Naf)