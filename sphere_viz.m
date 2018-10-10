% sphere visualizations
close all; 
N = 10;
[u,v] = meshgrid(linspace(0,2*pi,N),linspace(0,pi,N));
S = reshape([u,v],N*N,2);
P = [cos(S(:,1)).*sin(S(:,2)),sin(S(:,1)).*sin(S(:,2)),cos(S(:,2))]; 
% Jacobian of sphere parametrization
J = @(u,v) [-sin(u).*sin(v), cos(u).*cos(v); cos(u).*sin(v), sin(u).*cos(v); 0, -sin(v)];
% normalize points in domain
V = S./sqrt(sum(S.^2,2));
% or sample random directions
% V = 2*rand(N^2,2) - 1; V = V./sqrt(sum(V.^2,2));
% or a random constant direction
V = 2*rand(1,2) -1; V = repmat(V/norm(V),N^2,1);
% or a constant coordinate direction
% V = repmat([0,1],N^2,1);
% tangent vectors (pushforward)
Vt = zeros(N^2,3); for ii=1:N^2, Vt(ii,:) = V(ii,:)*J(S(ii,1),S(ii,2))'; end
Vt = Vt./sqrt(sum(Vt.^2,2));

% mesh the entire sphere
nmesh = 50; [XX,YY,ZZ] = sphere(nmesh);
fig = figure; 
subplot(1,2,1), scatter(S(:,1),S(:,2),'filled','k'); hold on; axis equal; grid on; axis([-1,2*pi+1,-1,pi+1]);
subplot(1,2,1), quiver(S(:,1),S(:,2),V(:,1),V(:,2),1,'k');
subplot(1,2,2), mesh(XX,YY,ZZ,ones(size(ZZ))); hold on; axis equal; colormap gray; alpha(0.75);
subplot(1,2,2), scatter3(P(:,1),P(:,2),P(:,3),'filled','k');
subplot(1,2,2), quiver3(P(:,1),P(:,2),P(:,3),Vt(:,1),Vt(:,2),Vt(:,3),1,'k');
fig.CurrentAxes.Visible = 'off';