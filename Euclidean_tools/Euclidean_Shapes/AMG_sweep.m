% Postprocess AMG sweep. Run SU2_postproc.m first!
clc; close all; clearvars;

addpath ~/AMG/Euclidean_tools/Euclidean_Shapes/
addpath ~/RESEARCH/SHDP/

datapath = '~/RESEARCH/SHDP/';
imgpath = '~/AMG/Euclidean_tools/Euclidean_Shapes/figs_AMG_flows/';

% load AMG_postproc
load([datapath,'AMG_postproc.mat']);

%% Read sweep data from sub-directory
cd(datapath)
[F] = readsu2_forces(dir('./forces'),3);

%% plot force results
set(0,'defaulttextInterpreter','latex')
fig = figure; plt = plot(t,F,'o-','linewidth',2); hold on; grid on;
% format
ylabel('$$C_d$$','fontsize',20)
xlabel('$$t$$, geodesic distance','fontsize',20)
set(gca,'FontSize',20)

%% determine warnings
dirstruct = dir([datapath,'WARNINGS']);
N_warn = size(dirstruct,1)-2; WARN = ones(N_warn,1);
for i=1:N_warn
    ind1 = strfind(dirstruct(i+2).name,'_');
    ind2 = strfind(dirstruct(i+2).name,'.');
    WARN(i) = str2double(dirstruct(i+2).name(ind1+1:ind2-1));
end
% plot warnings
scatter(t(WARN),F(WARN),200,'rx','linewidth',2)

%% plot in image version
imgfig = figure; dirstruct = dir(imgpath);
gifname = 'AMG_sweep.gif';
for i=1:length(t)
    figure(imgfig);
    hi = imshow([imgpath,dirstruct(i+2).name]); hold on;
    haxes = axes('pos',[0.125 0.6 .35 .35]); set(gca, 'color', 'none');
    plt  = plot(t,F,'o-','linewidth',2); hold on; grid on;
    sctr = scatter(t(WARN),F(WARN),200,'rx','linewidth',2);
    mv_sctr = scatter(t(i),F(i),50,'filled','k'); axis tight;
    % format
    ylabel('$$C_d$$','fontsize',14)
    xlabel('$$t$$, geodesic distance','fontsize',14)
    set(gca,'FontSize',14)
    
    % build gif
    figure(imgfig); frame = getframe(imgfig); 
    [A,map] = rgb2ind(frame2im(frame),256);
    if i == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.2);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
    end
    
    delete([mv_sctr,haxes,hi,plt,sctr]);
end


%% background image version
% imgfig = figure;
% for i=1:length(t)
%     figure(imgfig);
%     plt  = plot(t,F,'o-','linewidth',2); hold on; grid on;
%     sctr = scatter(t(WARN),F(WARN),200,'rx','linewidth',2);
%     mv_sctr = scatter(t(i),F(i),50,'filled','k'); axis tight;
%     % format
%     ylabel('$$C_d$$','fontsize',14)
%     xlabel('$$t$$, geodesic distance','fontsize',14)
%     set(gca,'FontSize',14)
%     % read in and display image
%     I=imread(dirstruct(i+2).name); I  = I(end:-1:1,:,:);
%     hi = image(0.9*xlim,0.95*ylim,I);
%     uistack(hi,'bottom')
%     % build gif
%     figure(imgfig); frame = getframe(imgfig); 
%     [A,map] = rgb2ind(frame2im(frame),256);
%     if i == 1
%         imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.2);
%     else
%         imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
%     end
%     
%     delete([mv_sctr,haxes,HIIMAGE,plt,sctr]);
% end
