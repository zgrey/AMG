% Compact embedded submanifold AMG (e.g., a sphere)
clc; close all; clearvars;

%% Convergence study
% convergence study values (NN and NT are amounts, N = 2.^NN and T = 2.^Tu are upper bounds)
NN = 10; NT = 1; Tu = 0; nboot = 100;
% combinations of N and T for convergence study
[Ni,Ti] = meshgrid(linspace(1,NN,NN),linspace(0,Tu,NT)); Ni = reshape(2.^Ni,NN*NT,1); Ti = reshape(2.^Ti,NN*NT,1);
% precondition metric vectors for convergence study
err = ones(NN*NT,nboot);

for j = 1:nboot
for i = 1:NN*NT
clc; fprintf('%0.2f%% complete...(%i / %i bootstrap)\n',i/(NN*NT)*100,j,nboot);
rng(j)
%% The sphere
% generate random values in parametrization domain
N = Ni(i);
% entire upper hemisphere
% r = 0.99*rand(N,1); th = 2*pi*rand(N,1);
% smaller neighborhoods
% r = (0.75-0.25)*rand(N,1)+ 0.25; th = pi/4*rand(N,1);
% r = (0.9-0.7)*rand(N,1)+ 0.7; th = pi/2*rand(N,1);
% r = 0.99*rand(N,1); th = pi/4*rand(N,1);
r = 0.25*rand(N,1); th = 2*pi*rand(N,1);

% sample random ball in parametrization domain
S = r.*[cos(th),sin(th)];
% sample random directions
V = 2*rand(N,2) - 1; V = V./sqrt(sum(V.^2,2));
% mesh the entire sphere
nmesh = 50; [XX,YY,ZZ] = sphere(nmesh);

% parametrize upper hemisphere
Zcurv = @(X,Y) sqrt(1-X.^2 - Y.^2); Z = Zcurv(S(:,1),S(:,2));
P = [S(:,1),S(:,2),Z];
% partials of sphere parametrization
Zx = @(X,Y) -X./(sqrt(1-X.^2 - Y.^2)); Zy = @(X,Y) -Y./(sqrt(1-X.^2 - Y.^2));
% Jacobian of sphere parametrization
J = @(X,Y) [repmat([1 0; 0 1],1,1,N); reshape(Zx(X,Y),1,1,N), reshape(Zy(X,Y),1,1,N)];
JJ = J(S(:,1),S(:,2));
% tangent vectors
Vt = zeros(N,3); for ii=1:N, Vt(ii,:) = V(ii,:)*JJ(:,:,ii)'; end
Vt = Vt./sqrt(sum(Vt.^2,2));
% exponential map
Exp = @(t,Vt,P) kron(P,cos(t)) + kron(Vt./sqrt(sum(Vt.^2,2)),sin(t));
% compute exponentials along random tangent vectors
k = 25; t = 0.5*linspace(-1,1,k)'; GG = reshape(Exp(t,Vt,P),k,N,3);
% logarithmic map
Log = @(p,P) acos(P*p').*(P - P*p'.*repmat(p,size(P,1),1))./sqrt(sum((P - P*p'.*repmat(p,size(P,1),1)).^2,2));

%% Ambient map on the sphere
m = 3; XY = [reshape(XX,(nmesh+1)^2,1),reshape(YY,(nmesh+1)^2,1),reshape(ZZ,(nmesh+1)^2,1)]; rng(47);

% linear ambient function
% a = 2*rand(m,1)-1; a = a/norm(a); Func = @(XY) XY*a; Grad = @(XY) repmat(a',size(XY,1),1);
% quadratic ambient ridge of rank(H) = r <= floor(m/2)
% r = 1; H = zeros(m); H(floor(m/2):floor(m/2)+r-1,floor(m/2):floor(m/2)+r-1) = eye(r); Func = @(X) sum((X*H).*X,2); Grad = @(X) X*H;
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

%% Karcher mean & Visualization
% initialize gradient descent optimization routine
p0  = Exp(0.01,Vt(1,:),P(1,:)); Jp = [1 0; 0 1; Zx(p0(1),p0(2)) Zy(p0(1),p0(2))]; vt = rand(1,2)*Jp'; vt = vt/norm(vt); v = 1;
% Fletcher et al. (Algorithm 1)
while norm(v) > 1e-6
    v = 1/N*sum(Log(p0,P),1);
    p0 = Exp(norm(v),v,p0); 
    Jp = [1 0; 0 1; Zx(p0(1),p0(2)) Zy(p0(1),p0(2))]; vt = rand(1,2)*Jp'; vt = vt/norm(vt);
end

%% Active Manifold-Geodesics
T = Ti(i);
% compute basis for tangent space (an orthogonal tangent vector)
vt2 = cross(p0,vt); Ptan = [p0 + vt; p0 + vt2; p0 - vt; p0 - vt2; p0 + vt;];

% compute AMG geodesic point set
k = 2; t = T*linspace(-1,1,k)'; Vg = sqrt(sum(Gt.^2,2)); Gset = zeros(k,N,3);
for ii=1:N, Gset(:,ii,:) = Exp(t*Vg(ii),Gt(ii,:),P(ii,:)); end
% refined geodesic samples (for visualization)
kk = 50; tt = T*linspace(-1,1,kk)'; GGset = zeros(kk,N,3);
for ii=1:N, GGset(:,ii,:) = Exp(tt*Vg(ii),Gt(ii,:),P(ii,:)); end

% logarithmic map of geodesic point set
Pset = reshape(Gset,k*N,3);
% logarithmic discrepancies of geodesic point set
Vlog1 = Log(p0,Pset(1:2:end-1,:));
Vlog2 = Log(p0,Pset(2:2:end,:));
Vlog = Vlog1 - Vlog2;
% SVD of tangential coordinates for logarithmic map of geodesic point set
[U,D,~] = svd(1/sqrt(N*(k-1)*(T^2+1))*[vt;vt2]*(Vlog)',0); U = U'*[vt;vt2];
% Remove bias from PGA
% logarithmic map of sampled point set
Vlogx = Log(p0,P);
[Ux,Dx,~] = svd(1/sqrt(N)*[vt;vt2]*Vlogx',0); Ux = Ux'*[vt;vt2];
% SVD of tangential coordinates for logarithmic map of geodesic point set
% [U,D,~] = svd(1/sqrt(N*(k-1)*(T^2+1))*[vt;vt2]*(Vlog-repmat(Vlogx,2,1))',0); U = U'*[vt;vt2];

% compute error w.r.t Mukherjee embedding definition
[Uemb,~,~] = svd(Grnd'); W = (eye(3) - p0'*p0)*Uemb(:,1); W = W./norm(W);
% subspace distance to the embedding definition
err(i,j) = norm(W*W' - U(1,:)'*U(1,:),2);

% compute active and inactive manifold-geodesic
AMG = Exp(2*tt,U(1,:),p0); IAMG = Exp(2*tt,U(2,:),p0);
% compute Mukherjee embedding definition geodesic
EmbG = Exp(2*tt,W',p0);

% compute AMG shadow plot
% project logarithmic map of random points
Gy = Log(p0,P)*U';
% first-singular projected geodesic points
Py = Exp(Gy(:,1),U(1,:),p0);

%% AMG level-set optimizatoin:
% find the point on the AMG that minimizes the function's variability over the IAMG
options = optimset('TolX',eps,'TolFun',eps,'Display','on');
obj = @(t) range(  Func( Exp(linspace(-10,10,500)',U(2,:),Exp(t,U(1,:),p0)) )  );
[tAMG,obj_opt,exitflag] = fminsearch(obj,0,options);
pAMG = Exp(tAMG,U(1,:),p0);
% compute inactive active manifold-geodesic at tAMG
IAMG_level = Exp(2*tt,U(2,:),pAMG);

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
end
end

if NT ~= 1
% compute subspace distance convergence rate
M = [ones(NT*NN,1) log10(Ti) log10(Ni)]; cerr = M \ log10(mean(err,2));
% coefficient of determination for convergence estimate
Rsq = 1 - sum((log10(mean(err,2)) - M*cerr).^2)/sum((log10(mean(err,2)) - mean(log10(mean(err,2)))).^2);
fprintf('Sub. dist. T convergence rate 10^%f (R^2 = %f)\n',cerr(2),Rsq);
fprintf('Sub. dist. N convergence rate 10^%f (R^2 = %f)\n',cerr(3),Rsq);
elseif NT == 1
% compute subspace distance convergence rate
M = [ones(NN,1) log10(Ni)]; cerr = M \ log10(mean(err,2));
% coefficient of determination for convergence estimate
Rsq = 1 - sum((log10(mean(err,2)) - M*cerr).^2)/sum((log10(mean(err,2)) - mean(log10(mean(err,2)))).^2);
fprintf('Sub. dist. N convergence rate 10^%f (R^2 = %f)\n',cerr(2),Rsq);
end

%% Visualizations
%% visualize function on sphere
fig = figure;
% simple filled surface
surf(XX,YY,ZZ,reshape(F,(nmesh+1),(nmesh+1)),'FaceColor','interp','FaceLighting','gouraud','EdgeAlpha',0);
% fancy filled surface
% surf(XX,YY,ZZ,reshape(F,(nmesh+1),(nmesh+1)),'FaceColor','interp','FaceLighting','gouraud','EdgeAlpha',0.25); alpha(0.55);
% contour style mesh
% mesh(XX,YY,ZZ,reshape(F,(nmesh+1),(nmesh+1)),'EdgeColor','interp'); alpha(0.75);
% no map, sphere only
% mesh(XX,YY,ZZ,ones(size(ZZ))); hold on; alpha(0); axis equal; colormap gray;
color = caxis; hold on; axis equal; colorbar;
% random points
scatter3(P(:,1),P(:,2),P(:,3),50,'filled','cdata',Frnd,'MarkerEdgeColor','k','linewidth',1);
% random tangent vectors
quiver3(P(:,1),P(:,2),P(:,3),Vt(:,1),Vt(:,2),Vt(:,3),1,'k','linewidth',2)
% exponential maps along random tangent vectors
for i=1:N, plot3(GG(:,i,1),GG(:,i,2),GG(:,i,3),'r','linewidth',2), end
% gradients of function at random points
% quiver3(P(:,1),P(:,2),P(:,3),Gt(:,1),Gt(:,2),Gt(:,3),1,'k','linewidth',2)
fig.CurrentAxes.Visible = 'off';

%% visualize converged Karcher mean
scatter3(p0(1),p0(2),p0(3),75,'k','filled','MarkerEdgeColor','k','linewidth',2);

%% visualize tangent space at mean
% quiver3(p0(1),p0(2),p0(3),vt(1),vt(2),vt(3),1,'k--','linewidth',2);
% quiver3(p0(1),p0(2),p0(3),vt2(1),vt2(2),vt2(3),1,'k--','linewidth',2);
% tangent plane at mean
plot3(Ptan(:,1),Ptan(:,2),Ptan(:,3),'k','linewidth',2)

%% plot geodesic set paths
% for i=1:N, plot3(GGset(:,i,1),GGset(:,i,2),GGset(:,i,3),'k','linewidth',1.5), end
% plot the geodesic set points
% scatter3(reshape(Gset(:,:,1),k*N,1),reshape(Gset(:,:,2),k*N,1),reshape(Gset(:,:,3),k*N,1),'k')

%% visualize log map vectors
% quiver3(repmat(p0(1),k*N,1),repmat(p0(2),k*N,1),repmat(p0(3),k*N,1),Vlog(:,1),Vlog(:,2),Vlog(:,3),1,'b','linewidth',1);
% visualize new active manifold-geodesic basis
quiver3(repmat(p0(1),2,1),repmat(p0(2),2,1),repmat(p0(3),2,1),U(:,1),U(:,2),U(:,3),1,'k','linewidth',2);
% visualize Mukherjee embedding-direction 
quiver3(p0(1),p0(2),p0(3),W(1),W(2),W(3),1,'r','linewidth',2);
% visualize PGA basis
quiver3(repmat(p0(1),2,1),repmat(p0(2),2,1),repmat(p0(3),2,1),Ux(:,1),Ux(:,2),Ux(:,3),1,'c','linewidth',2);

%% visualize active and inactive manifold-geodesic
plot3(AMG(:,1),AMG(:,2),AMG(:,3),'k','linewidth',2);
plot3(IAMG(:,1),IAMG(:,2),IAMG(:,3),'k--','linewidth',2);
plot3(EmbG(:,1),EmbG(:,2),EmbG(:,3),'r','linewidth',2);

%% visualize inactive active manifold-geodesic at tAMG
% scatter3(pAMG(1),pAMG(2),pAMG(3),75,'filled','r');
% plot3(IAMG_level(:,1),IAMG_level(:,2),IAMG_level(:,3),'r--','linewidth',2);

%% visualize rotating geodesics and AMG shadow
% line color scaling
Grscl = 10;
% for i=1:size(R,1)/Nr*Nr, plot3(Gr(:,i,1),Gr(:,i,2),Gr(:,i,3),'linewidth',1,'color',abs(1-Grscl^thr(i)/Grscl)*ones(3,1)), end
% geodesic sweeps
fig2 = figure; hold on;
for ii=1:size(R,1)/Nr*Nr, plot(tgr,Func(reshape(Gr(:,ii,:),Ngr,3)),'-','Color',abs(1-Grscl^thr(ii)/Grscl)*ones(3,1)); end
% replot AMG sweep
plot(tgr,Func(reshape(Gr(:,1,:),Ngr,3)),'k','linewidth',2);
% replot IAMG sweep
plot(tgr,Func(reshape(Gr(:,Nr,:),Ngr,3)),'k--','linewidth',2);
% scatter plot over inner products of Log projection
scatter(Gy(:,1),Frnd,50,'filled','cdata',Frnd); caxis(color);
ylabel 'f(Exp_x(t; v))'; xlabel 't'
% plot AMG level-set optimizaton inactive geodesic response
plot(tgr,Func(Exp(2*tgr,U(2,:),pAMG)),'r--','linewidth',2);

%% Convergence of subspace distance
if NN ~= 1 && NT ~= 1
% contour plot of linear fit to log10(err) for convergence rate estimates
figure; subplot(1,2,1), contourf(reshape(log10(Ni),NT,NN),reshape(log10(Ti),NT,NN),reshape([ones(NT*NN,1), log10(Ti), log10(Ni)]*cerr,NT,NN),15); hold on; axis square;
colorbar;
% colorbar('Ticks',[-16,-14,-12,-10,-8,-6,-4,-2,0]); caxis([-16,-1]); 
xlabel('$$\log_{10}(N)$$','Interpreter','latex'); ylabel('$$\log_{10}(T)$$','Interpreter','latex');
title(['$$\log_{10}\left(\Vert\hat{W_',num2str(r),'}\hat{W_',num2str(r),'}^T - \hat{U_',num2str(r),'}\hat{U_',num2str(r),'}^T\Vert_2\right)$$'],'Interpreter','latex');
% surface plot of log10(err) raw data
subplot(1,2,2), surf(reshape(log10(Ni),NT,NN),reshape(log10(Ti),NT,NN),reshape(log10(mean(err,2)),NT,NN));
colorbar; caxis([min(log10(mean(err,2))),max(log10(mean(err,2)))]);
hold on; shading interp; view([0,0,1]); axis([1,log10(max(Ni)),0,log10(max(Ti))]);
% colorbar('Ticks',[-16,-14,-12,-10,-8,-6,-4,-2,0]); caxis([-16,-1]); 
% contour expected rates over raw data
subplot(1,2,2), contour(reshape(log10(Ni),NT,NN),reshape(log10(Ti),NT,NN),reshape([ones(NT*NN,1), log10(Ti), log10(Ni)]*[cerr(1);2;0.5],NT,NN),15,'--','linecolor',0.5*ones(3,1));
xlabel('$$\log_{10}(N)$$','Interpreter','latex'); ylabel('$$\log_{10}(T)$$','Interpreter','latex'); axis square; axis([min(log10(Ni)),max(log10(Ni)),min(log10(Ti)),max(log10(Ti))]);
title(['$$\log_{10}\left(\Vert\hat{W_',num2str(r),'}\hat{W_',num2str(r),'}^T - \hat{U_',num2str(r),'}\hat{U_',num2str(r),'}^T\Vert_2\right)$$'],'Interpreter','latex');
elseif NN ~= 1
    figure; loglog(Ni,mean(err,2),'ko-','linewidth',2,'MarkerSize',8); hold on; grid on;
    err_lb = min(err,[],2); err_ub = max(err,[],2);
    h = fill([Ni; Ni(end:-1:1)],[err_lb; err_ub(end:-1:1)],0.5*ones(1,3)); h.FaceAlpha = 0.25;
    xlabel 'N'; ylabel 'subspace distance'
end
