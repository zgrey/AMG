% ralphie example
clc; close all; clearvars;
skp = 4; NN = 1000; G = [0 -1;1 0];

%% image processing
Rim = imread('./ralphie.png'); BW = imbinarize(Rim(:,:,1));
Rbnd = bwboundaries(~BW); bnd = Rbnd{3}*G; % use Rbnd{2}-inner or Rbnd{3}-outer
% select subset (ensure point set is from boundary of a simply-connected domain)
bndsub = bnd(1:skp:end,:); N = size(bndsub,1);
% shift and scale
bndsub = 2*(bndsub - min(bndsub,[],1))./(max(bndsub,[],1)-min(bndsub,[],1)) - 1;
x0 = [mean(bndsub(:,1)),mean(bndsub(:,2))]; bndsub = bndsub - x0;

%% circular embedding
s0 = atan2(bndsub(:,2),bndsub(:,1)); angi = find(s0 < 0);
s0(angi) = s0(angi) + 2*pi; s0 = s0/(2*pi);
a0 = bndsub(:,1).*cos(2*pi*s0) + bndsub(:,2).*sin(2*pi*s0);
% take only unique points
[s0,ia] = unique(s0,'stable');
a0 = a0(ia); Nu = length(s0);

%% 3D circular embedding
[PP3,n3] = embrep3(bndsub,NN,'uni'); 

%% increasing nearest-neighbor sort
lambda = 2; CW = 0;
[ind,S,n,S_rc,pp] = radial_emb_sort(s0,a0,NN,lambda,CW,1,1);
S = -S;
%% Visualize shape
fig = figure; hold on; axis equal;
plot([bndsub(:,1);bndsub(1,1)],[bndsub(:,2);bndsub(1,2)],'k','linewidth',4,'Color',[0.929,0.694,0.125]);
quiver(zeros(length(ind),1),zeros(length(ind),1),S(ind,1),S(ind,2),0,...
    'showarrowhead','off','Color',[0.502,0.502,0.502])
scatter(bndsub(:,1),bndsub(:,2),'k');
scatter(0,0,'filled','ko')
scatter(S(ind,1),S(ind,2),'r.')
axis equal; fig.CurrentAxes.Visible = 'off';
% 3D circular embedding
plot(PP3(:,1),PP3(:,2),'b--');

%% GIF
aaa = ppval(pp,linspace(0,1,NN)'); t  = linspace(0,max(aaa),100);
filename = 'ralphie.gif';
for k=1:length(t)
    
    % CCW
    PPP = [min([t(k)*ones(NN,1),aaa],[],2).*cos(linspace(0,1,NN)*2*pi)',...
           min([t(k)*ones(NN,1),aaa],[],2).*sin(linspace(0,1,NN)*2*pi)'];
    % CW
%     PPP = [min([t(k)*ones(NN,1),aaa],[],2).*cos(linspace(1,0,NN)*2*pi)',...
%            min([t(k)*ones(NN,1),aaa],[],2).*sin(linspace(1,0,NN)*2*pi)'];

    if k > 1
        delete(h);
    end
    h = plot(PPP(:,1),PPP(:,2),'k--','linewidth',2);
    
    frame = getframe(fig);
    [A,map] = rgb2ind(frame2im(frame),256);
    if k == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    
end

% plot(NEAR(:,1),NEAR(:,2),'k--')
%% Visualize circ.-embedding
figure; grid on; hold on; axis equal;
scatter(s0,a0,'k'); scatter(s0(ind),a0(ind),'r.');
plot(linspace(0,1,NN),ppval(pp,linspace(0,1,NN)),'k--');