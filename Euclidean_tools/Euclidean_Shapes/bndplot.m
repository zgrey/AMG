% plot a boundary function

function [h] = bndplot(fig,x,y,f)

figure(fig); hold on;
% cmap = colormap;
% con_min = min(f);
% con_max = max(f);
% ind_c = round((size(cmap,1)-1)*f/(con_max-con_min))+1;
% 
% c = round(1+(size(cmap,1)-1)*(f - min(f))/(max(f)-min(f)));
% 
% set(gca,'ColorOrder',cmap(ind_c));

% plot line segments on existing fig
% for k = 1:(length(x)-1)
%     line(x(k:k+1),y(k:k+1),f(k:k+1),'color',cmap(c(k),:),'linewidth',5)
% end
% view(2); 
% axis tight; colorbar; caxis([con_min con_max]);

%%
% n = length(f);
% p = plot(x,y,'r', 'LineWidth',5);
% 
% % modified jet-colormap
% cd = [uint8(jet(n)*255) uint8(f)].';
% 
% drawnow
% set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)

%% Prepare matrix data
xx=[x x]; yy=[y y]; zz=zeros(size(xx)); cc =[f f];

% draw the surface
h=surf(xx,yy,zz,cc,'EdgeColor','interp','FaceColor','none','linewidth',8) ;
colormap(parula); shading interp; view(2); colorbar; axis tight;