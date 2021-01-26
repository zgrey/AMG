% Grassmannian random gedesic wing
% Zach Grey
clc; close all; clearvars;
rng(42);

% addpath ~/AMG/Euclidean_tools/Euclidean_Shapes/
% addpath ~/RESEARCH/SHDP/
% addpath D:\AMG_DATA\GitHub_Backup\Euclidean_tools\Euclidean_Shapes
% load D:\AMG_DATA\GitHub_Backup\Euclidean_tools\Euclidean_Shapes\AMG_postproc.mat
% load /media/zgrey/46AFC5285FA7ADF9/AMG_DATA/Cd_adj/AMG_postproc.mat
% load D:\AMG_DATA\Cd_adj\AMG_postproc.mat
load D:\AMG_DATA\PGA_samples\BACKUP\qiqi_PGA_meshes.mat; P = Gr_Pts; Naf = 1000;

Nsec = 4;
alpha = [0;0;0;0];
scl = [1;0.75;0.5;0.25];

%% pre-process distances
subsmpl = randperm(Naf,Nsec)';

Pgeo = P(:,:,subsmpl);

% compute distance
d = zeros(Nsec,1);
for i = 2:Nsec
    d(i) = dGr_np(Pgeo(:,:,i-1),Pgeo(:,:,i));
end

%% build geodesics
L = linspace(0,sum(d),50*Nsec);
cumd = cumsum(d);
Nc = histc(L,cumd);
% rescale to wing length
cumd = 4*cumd;

fig = figure;
hold on; axis equal; P = []; Shape = zeros(Naf,2); len = zeros(Naf,1);
for i=1:Nsec-1
    clc; disp([num2str(i),'/',num2str(Nsec),' walks complete']);
    t = linspace(0,1,Nc(i));
       
    % plot nominal shapes
%     P0geo = Pgeo(:,:,i)*Minv_avg';
%     h0 = plot3(P0geo(:,1),cumd(i)*ones(Naf,1),P0geo(:,2),'linewidth',1); hold on;
    
    % compute direction
    [H] = Gr_log(Pgeo(:,:,i),Pgeo(:,:,i+1));
    for ii=1:length(t)
        % shape geodesic
        Gr_geo = Gr_exp(t(ii),Pgeo(:,:,i),H);
        Shape = Gr_geo*Minv_avg';
        len = ((1-t(ii))*cumd(i) + t(ii)*cumd(i+1))*ones(Naf,1);
        % visualize
        plot3(Shape(:,1),len,Shape(:,2),'k'); hold on; axis equal;
        fig.CurrentAxes.Visible = 'off';
    end
    
end