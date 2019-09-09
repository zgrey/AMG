% Grassmannian random gedesic walk
% Zach Grey
clc; close all; clearvars;
rng(42);

% addpath ~/AMG/Euclidean_tools/Euclidean_Shapes/
% addpath ~/RESEARCH/SHDP/
addpath D:\AMG_DATA\GitHub_Backup\Euclidean_tools\Euclidean_Shapes
load AMG_postproc.mat
Nsmpl = 10;

%% pre-process distances
subsmpl = randperm(Naf,Nsmpl)';

Pgeo = P(:,:,subsmpl);
% repeat last shape for visual continuity
Pgeo = cat(3,Pgeo,Pgeo(:,:,1)); subsmpl = [subsmpl; subsmpl(1)];

d = zeros(Nsmpl+1,1);

for i = 2:Nsmpl+1
    % compute distance
    d(i) = dGr_np(Pgeo(:,:,i-1),Pgeo(:,:,i));
end

%% build geodesics
L = linspace(0,sum(d),20*Nsmpl);
cumd = cumsum(d);
Nc = histc(L,cumd);

fig = figure; gifname = 'Rnd_geo.gif';
subplot(3,1,1), hd  = plot(cumd,ones(length(cumd),1),'o-','linewidth',2,'MarkerSize',10);
hold on; axis tight;
for i=1:Nsmpl
    clc; disp([num2str(i),'/',num2str(Nsmpl),' walks complete']);
    t = linspace(0,1,Nc(i));
    
    P0geo = sqrt(mean(n) - 1)*Pgeo(:,:,i)*Minv(:,:,subsmpl(i))';
    subplot(3,1,3), h0 = plot(P0geo(:,1),P0geo(:,2),'g--'); hold on;
    P0geo = sqrt(mean(n) - 1)*Pgeo(:,:,i+1)*Minv(:,:,subsmpl(i+1))';
    subplot(3,1,3), h01 = plot(P0geo(:,1),P0geo(:,2),'r--'); hold on;
    subplot(3,1,2), h1 = plot(Pgeo(:,1,i),Pgeo(:,2,i),'g--');
    hold on;
    subplot(3,1,2), h2 = plot(Pgeo(:,1,i+1),Pgeo(:,2,i+1),'r--');
    % compute direction
    [H] = Gr_log(Pgeo(:,:,i),Pgeo(:,:,i+1));
    % build geodesic *.gif
    for ii=1:length(t)
        % shape geodesic
        Gr_geo = Gr_exp(t(ii),Pgeo(:,:,i),H);
        Sgeo = Gr_geo*Minv_avg';

        % visualize
        subplot(3,1,1), hL1 = plot([0, cumd(i) + t(ii)*d(i+1)],[1,1],'k','linewidth',2);
        subplot(3,1,1), hL2 = scatter(cumd(i) + t(ii)*d(i+1),1,'k','filled','linewidth',2);
        fig.CurrentAxes.Visible = 'off';
        subplot(3,1,2), hGr = plot(Gr_geo(:,1),Gr_geo(:,2),'k','linewidth',1); hold on; axis equal;
        fig.CurrentAxes.Visible = 'off';
        subplot(3,1,3), hS  = plot(Sgeo(:,1),Sgeo(:,2),'k','linewidth',2); hold on; axis equal;
        fig.CurrentAxes.Visible = 'off';
        
        % build gif
        figure(fig); frame = getframe(fig); 
        [A,map] = rgb2ind(frame2im(frame),256);
        if i == 1
            imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.2);
        else
            imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
        end

        delete([hL1,hL2,hGr,hS]);
    end
    delete([h0,h01,h1,h2]);
end