% Compact embedded submanifold AMG (e.g., a sphere)
clc; close all; clearvars;

%% The sphere
% mesh the entire sphere (for visualization)
nmesh = 50; [XX,YY,ZZ] = sphere(nmesh);
% parametrize upper hemisphere
Zcurv = @(X,Y) sqrt(1-X.^2 - Y.^2);
% partials of sphere parametrization
Zx = @(X,Y) -X./(sqrt(1-X.^2 - Y.^2)); Zy = @(X,Y) -Y./(sqrt(1-X.^2 - Y.^2));
% Jacobian of sphere parametrization
J = @(X,Y) [repmat([1 0; 0 1],1,1,length(X)); reshape(Zx(X,Y),1,1,length(X)), reshape(Zy(X,Y),1,1,length(X))];
% exponential map
Exp = @(t,Vt,P) kron(P,cos(t)) + kron(Vt./sqrt(sum(Vt.^2,2)),sin(t));
% logarithmic map
Log = @(p,P) acos(P*p').*(P - P*p'.*repmat(p,size(P,1),1))./sqrt(sum((P - P*p'.*repmat(p,size(P,1),1)).^2,2));

%% Parametrize a neighborhood
% entire upper hemisphere
rad =@(N) 0.99*rand(N,1); theta =@(N) 2*pi*rand(N,1);
% smaller neighborhoods
% rad =@(N) (0.75-0.25)*rand(N,1)+ 0.25; theta =@(N) pi/4*rand(N,1);
% rad =@(N) (0.9-0.7)*rand(N,1)+ 0.7; theta =@(N) pi/2*rand(N,1);
% rad =@(N) 0.99*rand(N,1); theta =@(N) pi/4*rand(N,1);
% rad =@(N) 0.25*rand(N,1); theta =@(N) 2*pi*rand(N,1);

%% Ambient map on the sphere
m = 3; XY = [reshape(XX,(nmesh+1)^2,1),reshape(YY,(nmesh+1)^2,1),reshape(ZZ,(nmesh+1)^2,1)]; rng(47);
% linear ambient function
% a = 2*rand(m,1)-1; a = a/norm(a); Func = @(XY) XY*a; Grad = @(XY) repmat(a',size(XY,1),1);
% quadratic ambient ridge of rank(H) = r <= floor(m/2)
% r = 1; H = zeros(m); H(floor(m/2):floor(m/2)+r-1,floor(m/2):floor(m/2)+r-1) = eye(r); Func = @(X) sum((X*H).*X,2); Grad = @(X) 2*X*H;
% quadratic with preferential directions (my method is typically better)
H = diag(linspace(1,m,m)); Func = @(X) sum((X*H).*X,2); Grad = @(X) 2*X*H;
% highly nonlinear ridge
% a = 2*rand(m,1)-1; a = a/norm(a); Func = @(XY) sin(2*pi*XY*a) + cos(pi/2*XY*a); Grad = @(XY) kron(sum(pi*cos(pi*XY*a) - pi/2*sin(pi/2*XY*a),2),sum(a,2)');
% highly nonlinear approximate ridge
% w = 0.1; aa = 2*rand(m,1)-1; aa = aa/norm(aa); [A,~] = svd(aa); Func = @(XY) sum(sin(pi*XY*aa) + cos(pi/2*XY*aa),2) + w*sum((sin(pi*XY))*A(:,2:end),2); Grad = @(XY) kron(sum(pi*cos(pi*XY*aa) - pi/2*sin(pi/2*XY*aa),2),sum(aa,2)') + w*sum(pi*cos(pi*XY)*A(:,2:end),2)/(m-1); 
% highly nonlinear non-ridge
% aa = rand(m,1)-1; aa = aa/norm(aa); [A,~] = svd(aa); Func = @(XY) sum(sin(pi*XY*aa) + cos(pi/2*XY*aa),2) + sum((sin(pi*XY))*A(:,2:end),2); Grad = @(XY) kron(sum(pi*cos(pi*XY*aa) - pi/2*sin(pi/2*XY*aa),2),sum(aa,2)') + sum(pi*cos(pi*XY)*A(:,2:end),2)/(m-1);

%% Subpsace convergence study
% convergence study values (NN and NT are amounts, N = 2.^NN and T = 0.5.^NT are upper bounds)
NN = 7; NT = 1; nboot = 1;
% combinations of N and T for convergence study
[Ni,Ti] = meshgrid(linspace(1,NN,NN),linspace(1,NT,NT)); Ni = reshape(2.^Ni,NN*NT,1); Ti = reshape(0.1.^Ti,NN*NT,1);
% precondition metric vectors for convergence study
err = ones(NN*NT,nboot); err_emb = err;

%% Karcher mean & Visualization
N = max(Ni); r = rad(N); th = theta(N);
% sample random ball in parametrization domain
S = r.*[cos(th),sin(th)]; P = [S(:,1),S(:,2),Zcurv(S(:,1),S(:,2));];
% initialize gradient descent optimization routine
p0  = P(1,:); Jp = [1 0; 0 1; Zx(p0(1),p0(2)) Zy(p0(1),p0(2))]; vt = rand(1,2)*Jp'; vt = vt/norm(vt); v = 1;
% Fletcher et al. (Algorithm 1)
while norm(v) > 1e-8
    v = 1/(N-1)*sum(Log(p0,P(2:end,:)),1);
    p0 = Exp(norm(v),v,p0); 
    Jp = [1 0; 0 1; Zx(p0(1),p0(2)) Zy(p0(1),p0(2))]; vt = rand(1,2)*Jp'; vt = vt/norm(vt);
end

%% Run convergence study
for j = 1:nboot
for i = 1:NN*NT
%% print progress
clc; fprintf('%0.2f%% complete...(%i / %i bootstrap)\n',i/(NN*NT)*100,j,nboot);
rng(j)

%% generate random values in parametrization domain
N = Ni(i); r = rad(N); th = theta(N);
% sample random ball in parametrization domain
S = r.*[cos(th),sin(th)];
% sample random directions
V = 2*rand(N,2) - 1; V = V./sqrt(sum(V.^2,2));
% compute coordinates and pushforward
P = [S(:,1),S(:,2),Zcurv(S(:,1),S(:,2));];
JJ = J(S(:,1),S(:,2));
% tangent vectors
Vt = zeros(N,3); for ii=1:N, Vt(ii,:) = V(ii,:)*JJ(:,:,ii)'; end
Vt = Vt./sqrt(sum(Vt.^2,2));

% compute random function and ambient gradient evaluations
F = Func(XY); G =Grad(XY); Frnd = Func(P); Grnd = Grad(P);
% tangential gradient of ambient function
Gt = Grnd - bsxfun(@times,sum(Grnd.*P,2),P);

%% Active Manifold-Geodesics
T = Ti(i);
% compute basis for tangent space (an orthogonal tangent vector)
vt2 = cross(p0,vt); Ptan = [p0 + vt; p0 + vt2; p0 - vt; p0 - vt2; p0 + vt;];

% compute AMG geodesic point set
k = 2; t = T*linspace(-1,1,k)'; Vg = sqrt(sum(Gt.^2,2)); Gset = zeros(k,N,3);
for ii=1:N, Gset(:,ii,:) = Exp(t*Vg(ii),Gt(ii,:),P(ii,:)); end

% log-map of geodesic point set
Pset = reshape(Gset,k*N,3);

% compute coordinates for PGA basis
% logarithmic map of sampled point set
Vlogx = Log(p0,P);
[Ux,Dx,~] = svd(1/sqrt(N)*[vt;vt2]*Vlogx',0); Ux = Ux'*[vt;vt2];

% [method]
% [geodesic extensions] log-map of geodesic extensions (large T)
% Vlog = Log(p0,Pset);

% [central diff] log-map discrepancies of geodesic point set (ind. of T)
% Vlog1 = Log(p0,Pset(1:2:end-1,:));
% Vlog2 = Log(p0,Pset(2:2:end,:));
% Vlog = 1/(2*T)*(Vlog2 - Vlog1);

% [fwd diff] log-map discrepancis of geodesic points set (large T)
% Vlog1 = Log(p0,P);
% Vlog2 = Log(p0,Pset(N+1:end,:));
% Vlog = 1/T*(Vlog2 - Vlog1);

% [central diff-ladder] 
Vlog = zeros(N,3); Nrungs = 1e4;
for ii=1:N
    Vlog(ii,:) = diff_ladder(P(ii,:),p0,Gt(ii,:),Exp,Log,Nrungs);
end

% [Schild's ladder]
% Vlog = zeros(N,3); Nrungs = 100;
% for ii=1:N
%     Vlog(ii,:) = schilds_ladder(P(ii,:),p0,Gt(ii,:),Exp,Log,Nrungs);
% end

% [angle preserving isometry] NOT WORKING
% Vlog = zeros(N,3);
% for ii=1:N
%     Vlog(ii,:) = iso_angle(P(ii,:),p0,Gt(ii,:),Ux',Log);
% end

% SVD of tangential vectors
disp('Computing important directions...');
% compute in intrinsic dimension
% [U,D,~] = svd(1/sqrt(N)*Ux*Vlog',0); U = U'*Ux;
% compute in extrinsic dimension
[U,D,~] = svd(1/sqrt(N)*Vlog',0); U = U(:,1:2)';

% compute error w.r.t Mukherjee embedding definition
[Uemb,~,~] = svd(Grnd'); W = (eye(3) - p0'*p0)*Uemb(:,2); %W = W./norm(W);
if U(1,:)*W < 0, W = -W; end
% subspace distance to the embedding definition
if exist('a','var')
    Proj_a = (eye(3) - p0'*p0)*a; Proj_a = Proj_a/norm(Proj_a);
    err(i,j) = norm(U(1,:)'*U(1,:) - Proj_a*Proj_a',2);
else
    err(i,j) = norm(W*W' - U(1,:)'*U(1,:),2);
end

end
end

%% Compute convergence rates
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

%% Compute active and inactive manifold-geodesic
% refined geodesic samples (for visualization)
kk = 50; tt = linspace(-2,2,kk)'; GGset = zeros(kk,N,3);
for ii=1:N, GGset(:,ii,:) = Exp(tt*Vg(ii),Gt(ii,:),P(ii,:)); end
AMG = Exp(2*tt,U(1,:),p0); IAMG = Exp(2*tt,U(2,:),p0);
% compute Mukherjee embedding definition geodesic
EmbG = Exp(2*tt,W',p0);
% project logarithmic map of random points
Gy = Log(p0,P)*U';
% first-singular projected geodesic points
Py = Exp(Gy(:,1),U(1,:),p0);

%% AMG level-set search
% find the point on the AMG that minimizes the function's variability over the IAMG
options = optimset('TolX',eps,'TolFun',eps,'Display','on');
obj = @(t) range(  abs(Func( Exp(linspace(-2,2,5000)',U(2,:),Exp(t,U(1,:),p0)) ))  );
% optimization to determine the inactive submanifold
tAMG = zeros(100,1); obj_eval = 10*ones(100,1);
for ii=1:100
    [tAMG(i),obj_eval(i)] = fminbnd(obj,-2,2,options);
end
[~,opti] = min(obj_eval); tAMG = tAMG(opti);
% secant method to determine the inactive submanifold (for ridge function)
% tol = 0.1; t1 = 0; t2 = 1; iter = 0;
% while obj(t1) > tol && iter < 500, tn = t1 - obj(t1)*(t1 - t2)/( obj(t1) - obj(t2) ); t2 = t1; t1 = tn;  iter = iter + 1; end
% tAMG = t1;

% compute points over inactive geodesic
pAMG = Exp(tAMG,U(1,:),p0);
% compute inactive active manifold-geodesic at tAMG
IAMG_level = Exp(2*tt,U(2,:),pAMG);

%% Visualizations
% visualize function on sphere
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
ind = 1:N; if N > 500, ind = 1:500; end
scatter3(P(ind,1),P(ind,2),P(ind,3),50,'filled','cdata',Frnd(ind),'MarkerEdgeColor','k','linewidth',1);
% random gradients
quiver3(P(ind,1),P(ind,2),P(ind,3),Gt(ind,1),Gt(ind,2),Gt(ind,3),1,'k')
fig.CurrentAxes.Visible = 'off';

% visualize Karcher mean
scatter3(p0(1),p0(2),p0(3),75,'k','filled','MarkerEdgeColor','k','linewidth',2);
% visualize tangent space at mean
plot3(Ptan(:,1),Ptan(:,2),Ptan(:,3),'k','linewidth',2)

% plot geodesic set paths
% refined geodesic samples (for visualization)
% kk = 50; tt = T*linspace(-1,1,kk)'; GGset = zeros(kk,N,3);
% for ii=1:N, GGset(:,ii,:) = Exp(tt*Vg(ii),Gt(ii,:),P(ii,:)); end
% for i=1:N, plot3(GGset(:,i,1),GGset(:,i,2),GGset(:,i,3),'k','linewidth',1.5), end
% plot the geodesic set points
% scatter3(reshape(Gset(:,:,1),k*N,1),reshape(Gset(:,:,2),k*N,1),reshape(Gset(:,:,3),k*N,1),'k')

% visualize log map vectors
% quiver3(repmat(p0(1),N,1),repmat(p0(2),N,1),repmat(p0(3),N,1),Vlog(:,1),Vlog(:,2),Vlog(:,3),1,'b','linewidth',1);
% visualize new active manifold-geodesic basis
quiver3(p0(1),p0(2),p0(3),U(1,1),U(1,2),U(1,3),1,'k','linewidth',2);
quiver3(p0(1),p0(2),p0(3),U(2,1),U(2,2),U(2,3),1,'k--','linewidth',2);
% visualize Mukherjee embedding projection 
quiver3(p0(1),p0(2),p0(3),W(1),W(2),W(3),1,'r','linewidth',2);
quiver3(p0(1)*ones(3,1),p0(2)*ones(3,1),p0(3)*ones(3,1),Uemb(1,:)',Uemb(2,:)',Uemb(3,:)','r--','linewidth',2);
% visualize PGA basis
quiver3(repmat(p0(1),2,1),repmat(p0(2),2,1),repmat(p0(3),2,1),Ux(:,1),Ux(:,2),Ux(:,3),1,'color',0.5*ones(1,3),'linewidth',2);

% visualize active and inactive manifold-geodesic
plot3(AMG(:,1),AMG(:,2),AMG(:,3),'k','linewidth',2);
plot3(IAMG(:,1),IAMG(:,2),IAMG(:,3),'k--','linewidth',2);
plot3(EmbG(:,1),EmbG(:,2),EmbG(:,3),'r','linewidth',2);
% visualize inactive active manifold-geodesic at tAMG
plot3(IAMG_level(:,1),IAMG_level(:,2),IAMG_level(:,3),'--','linewidth',2,'color',0.5*ones(1,3));

% AMG shadow plots and visualize rotating geodesics
fig2 = figure; hold on;
% Rotating geodesics (verification of AMG and IAMG)
Nr = 50; tr = linspace(0,1,Nr)';
% geodesic sweep of first quadrant of the tangent space
R = kron(U(1,:),1-tr)  + kron(U(2,:),tr);
% geodesic sweep of two quadrants of the tangent space
% R = [R; kron(-U(1,:),tr)  + kron(U(2,:),1-tr)];
% normalize and construct inner-product color metric
R = R./sqrt(sum(R.^2,2)); thr = abs(sum(R.*repmat(U(1,:),size(R,1)/Nr*Nr,1),2));
% compute rotated geodesics
Ngr = 100; 
Tr = 2*max([abs(max(Gy(:,1))),abs(min(Gy(:,1)))]); 
tgr = linspace(-Tr,Tr,Ngr)';
% exponential maps along rotated tangent vectors
Gr = reshape(Exp(tgr,R,repmat(p0,size(R,1),1)),Ngr,size(R,1),3);
% line color scaling
Grscl = 5;
% geodesic sweeps
for ii=1:size(R,1), plot(tgr,Func(reshape(Gr(:,ii,:),Ngr,3)),'-','Color',abs(1-Grscl^thr(ii)/Grscl)*ones(3,1)); end
% replot AMG sweep
plot(tgr,Func(Exp(tgr,U(1,:),p0)),'k','linewidth',2);
% replot IAMG sweep at mean
plot(tgr,Func(Exp(tgr,U(2,:),p0)),'k--','linewidth',2);
% plot IAMG level set approximation
plot(tgr,Func(Exp(tgr,U(2,:),pAMG)),'--','linewidth',2,'color',0.5*ones(1,3));
% plot Muhkerjee embedding geodesic
plot(tgr,Func(Exp(tgr,W',p0)),'r','linewidth',2);
% scatter plot over inner products of Log projection (shadow plot)
scatter(Gy(:,1),Frnd,50,'filled','cdata',Frnd); caxis(color);
ylabel 'f(Exp_x(t; v))'; xlabel 't';

% Convergence of subspace distance
if NN ~= 1 && NT ~= 1
    % contour plot of linear fit to log10(err) for convergence rate estimates
    figure; subplot(1,2,1), contourf(reshape(log10(Ni),NT,NN),reshape(log10(Ti),NT,NN),reshape([ones(NT*NN,1), log10(Ti), log10(Ni)]*cerr,NT,NN),15); hold on; axis square;
    colorbar;
    % colorbar('Ticks',[-16,-14,-12,-10,-8,-6,-4,-2,0]); caxis([-16,-1]); 
    xlabel('$$\log_{10}(N)$$','Interpreter','latex'); ylabel('$$\log_{10}(T)$$','Interpreter','latex');
    title(['$$\log_{10}\left(\Vert\hat{W_',num2str(1),'}\hat{W_',num2str(1),'}^T - \hat{U_',num2str(1),'}\hat{U_',num2str(1),'}^T\Vert_2\right)$$'],'Interpreter','latex');
    % surface plot of log10(err) raw data
    subplot(1,2,2), surf(reshape(log10(Ni),NT,NN),reshape(log10(Ti),NT,NN),reshape(log10(mean(err,2)),NT,NN));
    colorbar; caxis([min(log10(mean(err,2))),max(log10(mean(err,2)))]);
    hold on; shading interp; view([0,0,1]); axis([1,log10(max(Ni)),log10(min(Ti)),log10(max(Ti))]);
    % colorbar('Ticks',[-16,-14,-12,-10,-8,-6,-4,-2,0]); caxis([-16,-1]); 
    % contour expected rates over raw data
    subplot(1,2,2), contour(reshape(log10(Ni),NT,NN),reshape(log10(Ti),NT,NN),reshape([ones(NT*NN,1), log10(Ti), log10(Ni)]*[cerr(1);0;0.5],NT,NN),15,'--','linecolor',0.5*ones(3,1));
    xlabel('$$\log_{10}(N)$$','Interpreter','latex'); ylabel('$$\log_{10}(T)$$','Interpreter','latex'); axis square; axis([min(log10(Ni)),max(log10(Ni)),min(log10(Ti)),max(log10(Ti))]);
    title(['$$\log_{10}\left(\Vert\hat{W_',num2str(1),'}\hat{W_',num2str(1),'}^T - \hat{U_',num2str(1),'}\hat{U_',num2str(1),'}^T\Vert_2\right)$$'],'Interpreter','latex');
elseif NN ~= 1
    figure; loglog(Ni,mean(err,2),'ko-','linewidth',2,'MarkerSize',8); hold on; grid on;
    err_lb = min(err,[],2); err_ub = max(err,[],2);
    h = fill([Ni; Ni(end:-1:1)],[err_lb; err_ub(end:-1:1)],0.5*ones(1,3)); h.FaceAlpha = 0.25;
    xlabel 'N'; ylabel 'subspace distance'
end
