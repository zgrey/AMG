% Grassmannian random gedesic walk
% Zach Grey
clc; close all; clearvars;
rng(42);

% addpath ~/AMG/Euclidean_tools/Euclidean_Shapes/
% addpath ~/RESEARCH/SHDP/
% addpath D:\AMG_DATA\GitHub_Backup\Euclidean_tools\Euclidean_Shapes
% load D:\AMG_DATA\GitHub_Backup\Euclidean_tools\Euclidean_Shapes\AMG_postproc.mat
load /media/zgrey/46AFC5285FA7ADF9/AMG_DATA/Cd_adj/AMG_postproc.mat

Nsmpl = 10;

%% pre-process distances
subsmpl = randperm(Naf,Nsmpl)'; Nsmpl = Nsmpl + 1; subsmpl = [subsmpl; 535];

Pgeo = P(:,:,subsmpl);
% repeat last shape for visual continuity
Pgeo = cat(3,Pgeo,Pgeo(:,:,1)); subsmpl = [subsmpl; subsmpl(1)];

d = zeros(Nsmpl+1,1);

for i = 2:Nsmpl+1
    % compute distance
    d(i) = dGr_np(Pgeo(:,:,i-1),Pgeo(:,:,i));
end

%% build geodesics
L = linspace(0,sum(d),50*Nsmpl);
cumd = cumsum(d);
Nc = histc(L,cumd);

fig = figure; gifname = 'Rnd_geo.gif';
% subplot(3,1,1), hd  = plot(cumd,ones(length(cumd),1),'o-','linewidth',2,'MarkerSize',10);
subplot(3,1,1), hd  = plot(cos(cumd/sum(d)*2*pi),sin(cumd/sum(d)*2*pi),'o','linewidth',2,'MarkerSize',10); hold on;
subplot(3,1,1), plot(cos(linspace(0,2*pi,100)),sin(linspace(0,2*pi,100)),'linewidth',2,'MarkerSize',10,'color',hd.Color);
hold on; axis tight; axis equal;
for i=1:Nsmpl
    clc; disp([num2str(i),'/',num2str(Nsmpl),' walks complete']);
    t = linspace(0,1,Nc(i));
    
    % plot accumulating distance over Grassmannian
%     subplot(3,1,1), hg = plot(cumd(i),1,'go','linewidth',2,'MarkerSize',10);
%     subplot(3,1,1), hs = plot(cumd(i+1),1,'ro','linewidth',2,'MarkerSize',10);
    subplot(3,1,1), hg = plot(cos(cumd(i)/sum(d)*2*pi),sin(cumd(i)/sum(d)*2*pi),'go','linewidth',2,'MarkerSize',10);
    subplot(3,1,1), hs = plot(cos(cumd(i+1)/sum(d)*2*pi),sin(cumd(i+1)/sum(d)*2*pi),'ro','linewidth',2,'MarkerSize',10);
    
    % plot nominal shapes
%     P0geo = sqrt(n(i) - 1)*Pgeo(:,:,i)*Minv(:,:,subsmpl(i))';
    P0geo = Pgeo(:,:,i)*Minv_avg';
    subplot(3,1,3), h0 = plot(P0geo(:,1),P0geo(:,2),'g','linewidth',2); hold on;
%     P0geo = sqrt(n(i) - 1)*Pgeo(:,:,i+1)*Minv(:,:,subsmpl(i+1))';
    P0geo = Pgeo(:,:,i+1)*Minv_avg';
    subplot(3,1,3), h01 = plot(P0geo(:,1),P0geo(:,2),'r','linewidth',2); hold on;
    
    % plot representative element on the Stiefel
    subplot(3,1,2), h1 = plot(Pgeo(:,1,i),Pgeo(:,2,i),'g','linewidth',2);
    hold on;
    subplot(3,1,2), h2 = plot(Pgeo(:,1,i+1),Pgeo(:,2,i+1),'r','linewidth',2);
    
    % compute direction
    [H] = Gr_log(Pgeo(:,:,i),Pgeo(:,:,i+1));
    % build geodesic *.gif
    for ii=1:length(t)
        % shape geodesic
        Gr_geo = Gr_exp(t(ii),Pgeo(:,:,i),H);
        Sgeo = Gr_geo*Minv_avg';

        % visualize
%         subplot(3,1,1), hL1 = plot([0, cumd(i) + t(ii)*d(i+1)],[1,1],'k','linewidth',2);
%         subplot(3,1,1), hL2 = scatter(cumd(i) + t(ii)*d(i+1),1,100,'k','filled','linewidth',2);
        tcirc = linspace(0,(cumd(i) + t(ii)*d(i+1))*2*pi/sum(d),100);
        subplot(3,1,1), hL1 = plot(cos(tcirc),sin(tcirc),'k','linewidth',2);
        tcirc = linspace(0,(cumd(i) + t(ii)*d(i+1))*2*pi/sum(d),100);
        subplot(3,1,1), hL2 = scatter(cos((cumd(i) + t(ii)*d(i+1))*2*pi/sum(d)),...
            sin((cumd(i) + t(ii)*d(i+1))*2*pi/sum(d)),100,'k','filled','linewidth',2);
        
        fig.CurrentAxes.Visible = 'off';
        subplot(3,1,2), hGr = plot(Gr_geo(:,1),Gr_geo(:,2),'k','linewidth',2); hold on; axis equal;
        fig.CurrentAxes.Visible = 'off';
        subplot(3,1,3), hS  = plot(Sgeo(:,1),Sgeo(:,2),'k','linewidth',2); hold on; axis equal;
        fig.CurrentAxes.Visible = 'off';
        
        % build gif
        figure(fig); frame = getframe(fig); 
        [A,map] = rgb2ind(frame2im(frame),256);
        if ii == 1 && i == 1
            imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.05);
        else
            imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.05);
        end

        delete([hL1,hL2,hGr,hS]);
    end
    delete([h0,h01,h1,h2,hg,hs]);
end