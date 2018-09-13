% post-process QiQi airfoil database to check radial convexity of shapes
clc; clearvars; close all;
% read airfoil
dirstruct  = dir('~/airfoil_coordinate_database/coordinates/');
Naf = size(dirstruct,1)-2;
ind_emb = ones(Naf,1); nonconvcount = ones(Naf,1);
for i=1:Naf
clc; fprintf('%0.2f%% complete...\n',i/(Naf)*100);
% read nominal points from file
P0 = importairfoil([dirstruct(2+i).folder,'/',dirstruct(2+i).name]);
% number of points
n = size(P0,1);
% sort for atan2 going LE -> TE -> LE (CCW begining at LE)
[~,LE_i] = min(P0(:,1));
P0 = [P0(LE_i:end,:);P0(1:LE_i-1,:)];

%% affine transformations (i.e., p  = M*p0   + a)
% shift and scale to [-1,1] box
pu = max(P0,[],1)'; pl = min(P0,[],1)';
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
Minv = inv(M); ainv = -Minv*a;  

% transform points
P = (P0*M + repmat(a',n,1));

% rotate to mean of first and last point
ang = atan2(P(1,2),P(1,1)); G = [cos(ang) -sin(ang); sin(ang) cos(ang)];
P = P*G;

%% embedding
% compute unique angles
s0 = atan2(P(:,2),P(:,1));
% scale to [0,1]
s0 = s0/(2*pi) + 1/2;

% circular
nu = [cos(2*pi*s0) , sin(2*pi*s0)];
% elliptical
% a = 4; b = 1; nu = [-a*sin(s0*2*pi), b*cos(s0*2*pi)]*[0 -1; 1 0];

a0 = P(:,1).*nu(:,1) + P(:,2).*nu(:,2);
% take only unique points
[s0,ia] = unique(s0,'stable');
a0 = a0(ia); Nu = length(s0);
% sort points based on embedding
N = 5000;
[ind,S,nu,S_rc,pp,s0,a0]= radial_emb_sort(s0,a0,N,10,0);
ind_emb(i) = sum(ind ~= 0);

%% convex check
k = convhull(P(:,1),P(:,2));
nonconvcount(i) = n - length(k);

%% visualize
subplot(1,2,1), h1 = scatter(P(:,1),P(:,2),25,'filled','k'); axis equal; hold on;
subplot(1,2,1), h2 = plot(P(:,1),P(:,2),'k',P(k,1),P(k,2),'r',S_rc(:,1),S_rc(:,2),'g','linewidth',2);
subplot(1,2,2), h3 = plot(s0,a0,'ko-',linspace(0,1,N),ppval(pp,linspace(0,1,N)),'g--'); grid on; hold on;
subplot(1,2,2), h4 = scatter(s0(ind(ind ~= 0)),a0(ind(ind ~= 0)),'filled','r');
h5 = annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Directory index i=',num2str(i),' ',dirstruct(2+i).name], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
print(['./figs/',num2str(i),'_',dirstruct(2+i).name,'.png'],'-dpng')
%% clear plot
delete([h1;h2;h3;h4;h5]);
end
fprintf('%f%% (%i / %i) are radially convex\n',(Naf - length(find(ind_emb)))/Naf*100,(Naf - length(find(ind_emb))),Naf)
fprintf('%f%% (%i / %i) are convex\n',(Naf - length(find(nonconvcount)))/Naf*100,(Naf - length(find(nonconvcount))),Naf)