% Compact embedded submanifold AMG (e.g., a sphere)
clc; close all; clearvars; rng(47);

%% The sphere
% generate random values in parametrization domain
N = 500;
% entire upper hemisphere
% r = 0.99*rand(N,1); th = 2*pi*rand(N,1);
% smaller neighborhood(s)
% r = (0.75-0.25)*rand(N,1)+ 0.25; th = pi/4*rand(N,1);
% r = (0.9-0.7)*rand(N,1)+ 0.7; th = pi/2*rand(N,1);
r = 0.99*rand(N,1); th = pi/4*rand(N,1);
% r = 0.25*rand(N,1); th = 2*pi*rand(N,1);

% sample random ball in parametrization domain
S = r.*[cos(th),sin(th)];
% sample random directions
V = 2*rand(N,2) - 1; V = V./sqrt(sum(V.^2,2));
% mesh the entire sphere
nmesh = 25; [XX,YY,ZZ] = sphere(nmesh);

% parametrize upper hemisphere
Zcurv = @(X,Y) sqrt(1-X.^2 - Y.^2); Z = Zcurv(S(:,1),S(:,2));
P = [S(:,1),S(:,2),Z];
% partials of sphere parametrization
Zx = @(X,Y) -X./(sqrt(1-X.^2 - Y.^2)); Zy = @(X,Y) -Y./(sqrt(1-X.^2 - Y.^2));
% Jacobian of sphere parametrization
J = @(X,Y) [repmat([1 0; 0 1],1,1,N); reshape(Zx(X,Y),1,1,N), reshape(Zy(X,Y),1,1,N)];
JJ = J(S(:,1),S(:,2));
% tangent vectors
Vt = zeros(N,3); for i=1:N, Vt(i,:) = V(i,:)*JJ(:,:,i)'; end
Vt = Vt./sqrt(sum(Vt.^2,2));
% exponential map
Exp = @(t,Vt,P) kron(P,cos(t)) + kron(Vt./sqrt(sum(Vt.^2,2)),sin(t));
% compute exponentials along random tangent vectors
k = 25; t = 0.5*linspace(-1,1,k)'; GG = reshape(Exp(t,Vt,P),k,N,3);
% logarithmic map
Log = @(p,P,vt) ( P - p )*([vt', cross(vt',p')]*( ([vt', cross(vt',p')]'*[vt', cross(vt',p')]) \ [vt', cross(vt',p')]' ))';

%% Ambient map on the sphere
m = 3; XY = [reshape(XX,(nmesh+1)^2,1),reshape(YY,(nmesh+1)^2,1),reshape(ZZ,(nmesh+1)^2,1)]; rng(47);

% highly nonlinear ridge
a = 2*rand(m,1)-1; a = a/norm(a); Func = @(XY) sin(2*pi*XY*a) + cos(pi/2*XY*a); Grad = @(XY) kron(sum(pi*cos(pi*XY*a) - pi/2*sin(pi/2*XY*a),2),sum(a,2)');
% highly nonlinear approximate ridge
% w = 0.1; a = 2*rand(m,1)-1; a = a/norm(a); [A,~] = svd(a); Func = @(XY) sum(sin(pi*XY*a) + cos(pi/2*XY*a),2) + w*sum((sin(pi*XY))*A(:,2:end),2); Grad = @(XY) kron(sum(pi*cos(pi*XY*a) - pi/2*sin(pi/2*XY*a),2),sum(a,2)') + w*sum(pi*cos(pi*XY)*A(:,2:end),2)/(m-1); 
% highly nonlinear non-ridge
% a = rand(m,1)-1; a = a/norm(a); [A,~] = svd(a); Func = @(XY) sum(sin(pi*XY*a) + cos(pi/2*XY*a),2) + sum((sin(pi*XY))*A(:,2:end),2); Grad = @(XY) kron(sum(pi*cos(pi*XY*a) - pi/2*sin(pi/2*XY*a),2),sum(a,2)') + sum(pi*cos(pi*XY)*A(:,2:end),2)/(m-1);

% Compute values
F = Func(XY); G =Grad(XY); Frnd = Func(P); Grnd = Grad(P);
% Tangential gradient of ambient function
Gt = Grnd - bsxfun(@times,sum(Grnd.*P,2),P);

%% Visualize function on sphere
fig = figure;
% filled surface
surf(XX,YY,ZZ,reshape(F,(nmesh+1),(nmesh+1)),'FaceColor','interp','FaceLighting','gouraud','EdgeAlpha',0.25); alpha(0.55);
% contour style mesh
% mesh(XX,YY,ZZ,reshape(F,(nmesh+1),(nmesh+1)),'EdgeColor','interp'); alpha(0.75);
% no map, sphere only
% mesh(XX,YY,ZZ,ones(size(ZZ))); hold on; alpha(0); axis equal; %colormap gray;
color = caxis; hold on; axis equal; colorbar;
% random points
scatter3(P(:,1),P(:,2),P(:,3),50,'filled','cdata',Frnd,'MarkerEdgeColor','k','linewidth',1);
% random tangent vectors
% quiver3(P(:,1),P(:,2),P(:,3),Vt(:,1),Vt(:,2),Vt(:,3),1,'k','linewidth',2)
% exponential maps along random tangent vectors
% for i=1:N, plot3(GG(:,i,1),GG(:,i,2),GG(:,i,3),'r','linewidth',2), end
% gradients of function at random points
quiver3(P(:,1),P(:,2),P(:,3),Gt(:,1),Gt(:,2),Gt(:,3),1,'k','linewidth',2)
fig.CurrentAxes.Visible = 'off';

%% Karcher mean & Visualization
% initialize gradient descent optimization routine
p0  = Exp(0,Vt(1,:),P(1,:)); Jp = [1 0; 0 1; Zx(p0(1),p0(2)) Zy(p0(1),p0(2))]; vt = rand(1,2)*Jp'; vt = vt/norm(vt); v = 1;
% Fletcher et al. (Algorithm 1)
while norm(v) > 1e-6
    v = 1/N*sum(Log(p0,P,vt),1);
    
%     scatter3(p0(1),p0(2),p0(3),35,'filled','b');
%     quiver3(p0(1),p0(2),p0(3),v(1),v(2),v(3),1,'b','linewidth',2);
    
    p0 = Exp(norm(v),v,p0); 
    Jp = [1 0; 0 1; Zx(p0(1),p0(2)) Zy(p0(1),p0(2))]; vt = rand(1,2)*Jp'; vt = vt/norm(vt);
      
%     pause(0.5)
end
% visualize converged Karcher mean
scatter3(p0(1),p0(2),p0(3),75,'r','filled','MarkerEdgeColor','k','linewidth',2);

%% Active Manifold-Geodesics
% known mean for upper hemisphere
% p0 = [0,0,1]; vt = [rand(1,2),0]; vt = vt/norm(vt);
% compute basis for tangent space (an orthogonal tangent vector)
vt2 = cross(p0,vt); Ptan = [p0 + vt; p0 + vt2; p0 - vt; p0 - vt2; p0 + vt;];
% visualize tangent space at mean
% quiver3(p0(1),p0(2),p0(3),vt(1),vt(2),vt(3),1,'k--','linewidth',2);
% quiver3(p0(1),p0(2),p0(3),vt2(1),vt2(2),vt2(3),1,'k--','linewidth',2);
% tangent plane at mean
plot3(Ptan(:,1),Ptan(:,2),Ptan(:,3),'k','linewidth',2)

% compute & visualize AMG geodesic point set
T = 1; k = 2; t = T*linspace(-1,1,k)'; Vg = sqrt(sum(Gt.^2,2)); Gset = zeros(k,N,3);
for i=1:N, Gset(:,i,:) = Exp(t*Vg(i),Gt(i,:),P(i,:)); end
% refined geodesic samples (for visualization)
kk = 25; tt = T*linspace(-1,1,kk)'; GGset = zeros(kk,N,3);
for i=1:N, GGset(:,i,:) = Exp(tt*Vg(i),Gt(i,:),P(i,:)); end
% plot geodesic set paths
% for i=1:N, plot3(GGset(:,i,1),GGset(:,i,2),GGset(:,i,3),'k','linewidth',1.5), end
% plot the geodesic set points
% scatter3(reshape(Gset(:,:,1),k*N,1),reshape(Gset(:,:,2),k*N,1),reshape(Gset(:,:,3),k*N,1),'k')

% logarithmic map of geodesic point set
Vlog = Log(p0,reshape(Gset,k*N,3),vt);
% SVD of tangential coordinates for logarithmic map of geodesic point set
[U,D,~] = svd(1/sqrt(N*(k-1)*(T^2+1))*[vt;vt2]*Vlog',0); U = U'*[vt;vt2];
% visualize log map vectors
quiver3(repmat(p0(1),k*N,1),repmat(p0(2),k*N,1),repmat(p0(3),k*N,1),Vlog(:,1),Vlog(:,2),Vlog(:,3),1,'b','linewidth',1);
% visualize new active manifold-geodesic basis
quiver3(repmat(p0(1),2,1),repmat(p0(2),2,1),repmat(p0(3),2,1),U(:,1),U(:,2),U(:,3),1,'r','linewidth',2);

% compute & visualize active manifold-geodesic
AMG = Exp(2*tt,U(1,:),p0); IAMG = Exp(2*tt,U(2,:),p0);
plot3(AMG(:,1),AMG(:,2),AMG(:,3),'r','linewidth',2);
plot3(IAMG(:,1),IAMG(:,2),IAMG(:,3),'r--','linewidth',2);

% compute & visualize AMG shadow plot
% project logarithmic map of random points
Gy = Log(p0,P,vt)*U';
% visualize first-singular projected geodesic points
Py = Exp(Gy(:,1),U(1,:),p0);
plot3(Py(:,1),Py(:,2),Py(:,3),'ro','linewidth',2);

%% Rotating geodesics (verification of AMG and IAMG)
Nr = 25; tr = linspace(0,1,Nr)';
% geodesic sweep of first quadrant of the tangent space
R = kron(U(1,:),1-tr)  + kron(U(2,:),tr);
% geodesic sweep of two quadrants of the tangent space
% R = [R; kron(-U(1,:),tr)  + kron(U(2,:),1-tr)];
% normalize and construct inner-product color metric
R = R./sqrt(sum(R.^2,2)); thr = abs(sum(R.*repmat(U(1,:),size(R,1)/Nr*Nr,1),2));
% plot rotated tangent vectors
% figure(fig); quiver3(repmat(p0(1),size(R,1)/Nr*Nr,1), repmat(p0(2),size(R,1)/Nr*Nr,1), repmat(p0(3),size(R,1)/Nr*Nr,1), R(:,1), R(:,2), R(:,3),'Color',0.75*ones(3,1))
% compute rotated geodesics
Ngr = 100; Tr = 2*max([abs(max(Gy(:,1))),abs(min(Gy(:,1)))]); tgr = linspace(-Tr,Tr,Ngr)';
% exponential maps along rotated tangent vectors
Gr = reshape(Exp(tgr,R,repmat(p0,size(R,1)/Nr*Nr,1)),Ngr,size(R,1)/Nr*Nr,3);
% line color scaling
Grscl = 10;
% exponential map plots
% for i=1:size(R,1)/Nr*Nr, plot3(Gr(:,i,1),Gr(:,i,2),Gr(:,i,3),'linewidth',1,'color',abs(1-Grscl^thr(i)/Grscl)*ones(3,1)), end
% geodesic sweeps
fig2 = figure; hold on;
for i=1:size(R,1)/Nr*Nr, plot(tgr,Func(reshape(Gr(:,i,:),Ngr,3)),'-','Color',abs(1-Grscl^thr(i)/Grscl)*ones(3,1)); end
% replot AMG sweep
plot(tgr,Func(reshape(Gr(:,1,:),Ngr,3)),'k','linewidth',1);
% replot IAMG sweep
plot(tgr,Func(reshape(Gr(:,Nr,:),Ngr,3)),'k--','linewidth',1);
% scatter plot over inner products of Log projection
scatter(Gy(:,1),Frnd,50,'filled','cdata',Frnd); caxis(color);
ylabel 'f(Exp_x(t; v))'; xlabel 't'