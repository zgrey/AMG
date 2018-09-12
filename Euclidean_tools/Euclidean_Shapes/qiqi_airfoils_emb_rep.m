% post-process QiQi airfoil database to check radial convexity of shapes
clc; clearvars; close all;
% read airfoil
dirstruct  = dir('~/airfoil_coordinate_database/coordinates/');
Naf = size(dirstruct,1)-2;
ind_circ = ones(Naf,1); ind_ellp = ones(Naf,1); nonconvcount = ones(Naf,1);
for i=34%1:Naf
% read nominal points from file
P0 = importairfoil([dirstruct(2+i).folder,'/',dirstruct(2+i).name]);
% number of points
n = size(P0,1);

%% shift, scale, and rotate
% shift and scale to [-1,1] box
pu = max(P0,[],1)'; pl = min(P0,[],1)';
M    = diag(2./(pu-pl)); a    = -(M*pl + ones(2,1)); % i.e., x  = M*x0   + a    in [-1,1]_m
Minv = diag((pu-pl)./2); ainv = -Minv*a;             % i.e., x0 = Minv*x + ainv in [xl,xu]_m
% shift, and scale
P = (P0*M' + repmat(a',n,1));
% rotate to mean of first and last point
ang = atan2(P(1,2)/2 + P(end,2)/2,P(1,1)/2 + P(end,1)/2); G = [cos(ang) -sin(ang); sin(ang) cos(ang)];
P = P*G;

%% circular embedding
s0 = atan2(P(:,2),P(:,1)); s0(1) = 0; angi = find(s0 < 0);
s0(angi) = s0(angi) + 2*pi; s0 = s0/(2*pi);
% create monotone distribution of points
s0 = cumsum(s0);
nu = [cos(2*pi*s0) , sin(2*pi*s0)];
a0 = P(:,1).*nu(:,1) + P(:,2).*nu(:,2);
% take only unique points
[s0,ia] = unique(s0,'stable');
a0 = a0(ia); Nu = length(s0);
% sort points based on embedding
ind = radial_emb_sort(s0,a0,100,2,0);
ind_circ(i) = sum(ind ~= 0);

%% elliptical embedding
% a = 4; b = 1;
% s0 = atan2(P(:,2),P(:,1)); angi = find(s0 < 0);
% s0(angi) = s0(angi) + 2*pi; s0 = s0/(2*pi);
% nu = [-a*sin(s0*2*pi), b*cos(s0*2*pi)]*[0 -1; 1 0];
% a0 = P(:,1).*nu(:,1) + P(:,2).*nu(:,2);
% % take only unique points
% [s0,ia] = unique(s0,'stable');
% a0 = a0(ia); Nu = length(s0);
% ind = radial_emb_sort(s0,a0,100,2,0);
% ind_ellp(i) = sum(ind ~= 0);

%% convex check
k = convhull(P(:,1),P(:,2));
nonconvcount(i) = n - length(k);

%% visualize
subplot(1,2,1), h1 = scatter(P(:,1),P(:,2),25,'filled','k'); axis equal; hold on;
subplot(1,2,1), h2 = plot(P(:,1),P(:,2),'k',P(k,1),P(k,2),'r','linewidth',2);
subplot(1,2,2), h3 = plot(s0,a0,'ko-');
pause(0.25);
%% clear plot
% delete([h1;h2;h3]);
end