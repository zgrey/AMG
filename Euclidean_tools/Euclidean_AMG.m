%% Euclidean AMG demo script & convergence study
% reference: Grey, Z.J., Constantine, P.G., "A Riemannian View on Active Subspaces", TBD
close all; clearvars;
rng(47);

%% Things to modify
% dimension of Euclidean space (m > 1)
m = 10;
% convergence study values (NN and NT are amounts, 2.^Nu and 2.^Tu are upper bounds)
NN = 10; NT = 100; Nu = 10; Tu = 10;

%% Convergence study
% combinations of N and T for convergence study
[N,T] = meshgrid(linspace(1,Nu,NN),linspace(1,Tu,NT)); N = reshape(2.^N,NN*NT,1); T = reshape(2.^T,NN*NT,1);
% precondition metric vectors for convergence study
ang1 = 90*ones(NN*NT,1); dot1 = zeros(NN*NT,1); err = ones(NN*NT,1); ang2 = ang1;

for i=1:NN*NT
clc; fprintf('%0.2f%% complete...\n',i/(NN*NT)*100);
% uniform random samples
X = -1 + 2*rand(N(i),m);

%% Various real scalar-valued maps and gradient, uncomment line of interest:
% oscillating ridge
a = 2*rand(m,1)-1; a = a/norm(a); F = sin(2*pi*X*a) + cos(pi/2*X*a); Grad = kron(sum(pi*cos(pi*X*a) - pi/2*sin(pi/2*X*a),2),sum(a,2)');
% oscillating approximate ridge
% a = 2*rand(m,1)-1; a = a/norm(a); [A,~] = svd(a); F = sum(sin(pi*X*a) + cos(pi/2*X*a),2) + 0.1*sum((sin(pi*X))*A(:,2:end),2); Grad = kron(sum(pi*cos(pi*X*a) - pi/2*sin(pi/2*X*a),2),sum(a,2)') + 0.1*sum(pi*cos(pi*X)*A(:,2:end),2)/(m-1);
% oscillating non-ridge
% a = 2*rand(m,1)-1; a = a/norm(a); [A,~] = svd(a); F = sum(sin(pi*X*a) + cos(pi/2*X*a),2) + sum((sin(pi*X))*A(:,2:end),2); Grad = kron(sum(pi*cos(pi*X*a) - pi/2*sin(pi/2*X*a),2),sum(a,2)') + sum(pi*cos(pi*X)*A(:,2:end),2)/(m-1);
% exponential ridge
% a = 2*rand(m,1)-1; F = exp(X*a); Grad = repmat(a',N(i),1).*F;
% quadratic non-ridge with increasing parameter importance and rank(C) = m-1
% H = diag(linspace(0,1,m)); F = sum((X*H).*X,2); Grad = X*H;
% linear-quadratic ridge of rank(C) =  r <= m
% r = m; a = 2*rand(m,1)-1; H = [eye(r),zeros(r,m-r);zeros(m-r,m)]; F = sum((X*H).*X + X*a,2); Grad = X*H + repmat(a',N(i),1);
% quadratic ridge of rank(C) = r <= floor(m/2)
% r = 2; H = zeros(m); H(floor(m/2):floor(m/2)+r-1,floor(m/2):floor(m/2)+r-1) = eye(r); F = sum((X*H).*X,2); Grad = X*H;

%% Euclidean Active Manifold-Geodesics approximation
% k=1 geodesic extensions
t = [-T(i), T(i)];
% collate geodesic point set
G = kron(t',Grad) + repmat(X,2,1);
% svd of geodesic point set
[U,Sig,~] = svd(1/sqrt(N(i)*(T(i)^2+1))*G',0); U1 = U(:,1); U2 = U(:,2);

%% Active Subspace approximation
[W,D,~] = svd(1/sqrt(N(i))*Grad',0); W1 = W(:,1); W2 = W(:,2);

%% Convergence metrics
% convergence metric: subspace angles
ang1(i) = acos(abs(W1'*U1))*180/pi; ang2(i) = acos(abs(W2'*U2))*180/pi;
% convergence metric: inner product
dot1(i) = 1-abs(W1'*U1);
% convergence metric: subspace distance
% determine largest eigenvalue gap
lambda = diag(D); [~,r] = max(abs(lambda(2:end)-lambda(1:end-1)));
% subspace distance
err(i) = norm(W(:,1:r)*W(:,1:r)' - U(:,1:r)*U(:,1:r)',2);

end

%% Report subspace angles
% print subspace angles for the first singular vectors
% first values obtained
fprintf('1st u_1 angle = %f deg.\n',ang1(1));
fprintf('1st u_2 angle = %f deg.\n',ang2(1));
% last values obtained
fprintf('Last u_1 angle = %f deg.\n',ang1(end));
fprintf('Last u_2 angle = %f deg.\n',ang2(end));

%% Compute convergence rates
% screen for truncation/rounding errors
% determine negative truncation/roundoff and zeros to fill in log10 contours
dot1(dot1 <= 0) = eps;
% compute subspace distance convergence rate
M = [ones(NT*NN,1) log10(T) log10(N)]; cerr = M \ log10(err); cerr0 = [cerr(1); 2; 0.5];
% coefficient of determination for convergence estimate
Rsq = 1 - sum((log10(err) - M*cerr).^2)/sum((log10(err) - mean(log10(err))).^2);
fprintf('Sub. dist. T convergence rate 10^%f (R^2 = %f)\n',cerr(2),Rsq);
fprintf('Sub. dist. N convergence rate 10^%f (R^2 = %f)\n',cerr(3),Rsq);

% compute first singular vector inner product convergennce rate
cdot1 = M \ log10(dot1); cdot10 = [cdot1(1); 2; 0.5];
% coefficient of determination for convergence estimate
Rsq_dot1 = 1 - sum((log10(dot1) - M*cdot1).^2)/sum((log10(dot1) - mean(log10(dot1))).^2);
fprintf('First eig. vec. T convergence rate 10^%f (R^2 = %f)\n',cdot1(2), Rsq_dot1);
fprintf('First eig. vec. N convergence rate 10^%f (R^2 = %f)\n',cdot1(3), Rsq_dot1);

%% Visualizations
%% Visualizations of function
% 2D domain visualizations
if NT*NN < 1500 && m==2
% visualize function, gradients, and samples
subplot(1,2,1), scatter(X(:,1),X(:,2),50,'filled','cdata',F);
hold on; axis equal; xlabel('$$x_1$$','Interpreter','latex'); ylabel('$$x_2$$','Interpreter','latex');
subplot(1,2,1), quiver(X(:,1),X(:,2),Grad(:,1),Grad(:,2),'k'); grid on;
% subspace angle title block
gcf; title(['$$\theta_1 = ',num2str(ang1(end,end)),'$$ deg.'],'Interpreter','latex');

% plot AMG basis
subplot(1,2,1), quiver(0,0,U1(1),U1(2),'r','linewidth',2); 
subplot(1,2,1), quiver(0,0,U2(1),U2(2),'r--','linewidth',2);
% plot active subspace basis
subplot(1,2,1), quiver3(0,0,0,W1(1),W1(2),W1(3),'k--','linewidth',2); 
subplot(1,2,1), quiver3(0,0,0,W2(1),W2(2),W2(3),'k--','linewidth',2);
% scatter AMG geodesic point set (function stencil)
subplot(1,2,2), scatter(G(:,1),G(:,2),'k.');
axis equal; fig = gcf; fig.CurrentAxes.Visible = 'off';

% 3D domain visualizations
elseif NT*NN < 1500 && m==3
% visualize function, gradients, and samples
subplot(1,2,1), scatter3(X(:,1),X(:,2),X(:,3),50,'filled','cdata',F);
hold on; axis equal; xlabel('$$x_1$$','Interpreter','latex'); ylabel('$$x_2$$','Interpreter','latex'); zlabel('$$x_3$$','Interpreter','latex');
subplot(1,2,1), quiver3(X(:,1),X(:,2),X(:,3),Grad(:,1),Grad(:,2),Grad(:,3),'k'); grid on;
% subspace angle title block
gcf; title(['$$\theta_1 = ',num2str(ang1(end,end)),'$$ deg.'],'Interpreter','latex');

% plot AMG basis
subplot(1,2,1), quiver3(0,0,0,U1(1),U1(2),U1(3),'r','linewidth',2); 
subplot(1,2,1), quiver3(0,0,0,U2(1),U2(2),U2(3),'r--','linewidth',2);
% plot active subspace basis
subplot(1,2,1), quiver3(0,0,0,W1(1),W1(2),W1(3),'k--','linewidth',2); 
subplot(1,2,1), quiver3(0,0,0,W2(1),W2(2),W2(3),'k--','linewidth',2);
% scatter AMG geodesic point set (function stencil)
subplot(1,2,2), scatter3(G(:,1),G(:,2),G(:,3),'k.'); axis equal; fig = gcf; fig.CurrentAxes.Visible = 'off';
end

%% Shadow plots
% first principal direction shadows
figure; subplot(1,2,1), scatter(X*U1,F,'filled');
title('Active Manifold-Geodesics','Interpreter','latex'); xlabel('$$\hat{u}_1^Tx$$','Interpreter','latex'); ylabel('f','Interpreter','latex')
subplot(1,2,2), scatter(X*W1,F,'filled');
title('Active Subspaces','Interpreter','latex'); xlabel('$$\hat{w}_1^Tx$$','Interpreter','latex'); ylabel('f','Interpreter','latex')
% first and second principal direction shadows
if m > 2
figure; subplot(1,2,1), scatter(X*U1,X*U2,'filled','cdata',F);
title('Active Manifold-Geodesics','Interpreter','latex'); axis equal; xlabel('$$\hat{u}_1^Tx$$','Interpreter','latex'); ylabel('$$\hat{u}_2^Tx$$','Interpreter','latex'); colorbar;
subplot(1,2,2), scatter(X*W1,X*W2,'filled','cdata',F);
title('Active Subspaces','Interpreter','latex'); axis equal; xlabel('$$\hat{w}_1^Tx$$','Interpreter','latex'); ylabel('$$\hat{w}_2^Tx$$','Interpreter','latex'); colorbar;
end

%% Eigenvalues and 1st principal eigenvectors
% eigenvalues
figure; subplot(1,2,1), semilogy(diag(Sig),'o-','MarkerSize',10); hold on; plot(diag(D),'k*-'); grid on;
xlabel('i','interpreter','latex'); ylabel('$$\hat{\sigma}^i \, (Circles),\,\, \hat{\lambda}^i\, (Asterisk)$$','interpreter','latex'); title('Eigenvalues','Interpreter','latex');
% eigenvectors
subplot(1,2,2), stem(U1,'MarkerSize',10); hold on; stem(W1,'k*'); grid on;
xlabel('i','interpreter','latex'); ylabel('$$\hat{u}^i_1 \, (Circles),\,\, \hat{w}^i_1\, (Asterisk)$$','interpreter','latex'); title('$$1^{st}$$ Singular Eigenvectors','Interpreter','latex');

%% Convergence of first singular vector
if length(N) ~= 1
% contour plot of linear fit to log10(dot1) for convergence rate estimates
figure; subplot(1,2,1), contourf(reshape(log10(N),NT,NN),reshape(log10(T),NT,NN),reshape([ones(NT*NN,1), log10(T), log10(N)]*cdot1,NT,NN),15);
colorbar('Ticks',[-16,-14,-12,-10,-8,-6,-4,-2,0]); caxis([-16,-1]); hold on; axis equal;
xlabel('$$\log_{10}(N)$$','Interpreter','latex'); ylabel('$$\log_{10}(T)$$','Interpreter','latex'); title('$$\log_{10}(1-|\hat{w}_1^T\hat{u}_1|)$$','Interpreter','latex');
% surface plot of log10(dot1) raw data
subplot(1,2,2), surf(reshape(log10(N),NT,NN),reshape(log10(T),NT,NN),reshape(log10(dot1),NT,NN));
hold on; colorbar('Ticks',[-16,-14,-12,-10,-8,-6,-4,-2,0]); caxis([-16,-1]); shading interp; view([0,0,1]); axis([1,log10(max(N)),1,log10(max(T))]);
% contour expected rates over raw data
subplot(1,2,2), contour(reshape(log10(N),NT,NN),reshape(log10(T),NT,NN),reshape([ones(NT*NN,1), log10(T), log10(N)]*[cdot1(1);2;0.5],NT,NN),15,'--','linecolor',0.25*ones(3,1));
xlabel('$$\log_{10}(N)$$','Interpreter','latex'); ylabel('$$\log_{10}(T)$$','Interpreter','latex'); axis square; axis([min(log10(N)),max(log10(N)),min(log10(T)),max(log10(T))]);
title('$$\log_{10}(1-|\hat{w}_1^T\hat{u}_1|)$$','Interpreter','latex');
end

%% Convergence of subspace distance
if length(N) ~= 1
% contour plot of linear fit to log10(err) for convergence rate estimates
figure; subplot(1,2,1), contourf(reshape(log10(N),NT,NN),reshape(log10(T),NT,NN),reshape([ones(NT*NN,1), log10(T), log10(N)]*cerr,NT,NN),15);
colorbar('Ticks',[-16,-14,-12,-10,-8,-6,-4,-2,0]); caxis([-16,-1]); hold on; axis square;
xlabel('$$\log_{10}(N)$$','Interpreter','latex'); ylabel('$$\log_{10}(T)$$','Interpreter','latex');
title(['$$\log_{10}\left(\Vert\hat{W_',num2str(r),'}\hat{W_',num2str(r),'}^T - \hat{U_',num2str(r),'}\hat{U_',num2str(r),'}^T\Vert_2\right)$$'],'Interpreter','latex');
% surface plot of log10(err) raw data
subplot(1,2,2), surf(reshape(log10(N),NT,NN),reshape(log10(T),NT,NN),reshape(log10(err),NT,NN));
hold on; colorbar('Ticks',[-16,-14,-12,-10,-8,-6,-4,-2,0]); caxis([-16,-1]); shading interp; view([0,0,1]); axis([1,log10(max(N)),1,log10(max(T))]);
% contour expected rates over raw data
subplot(1,2,2), contour(reshape(log10(N),NT,NN),reshape(log10(T),NT,NN),reshape([ones(NT*NN,1), log10(T), log10(N)]*[cerr(1);2;0.5],NT,NN),15,'--','linecolor',0.5*ones(3,1));
xlabel('$$\log_{10}(N)$$','Interpreter','latex'); ylabel('$$\log_{10}(T)$$','Interpreter','latex'); axis square; axis([min(log10(N)),max(log10(N)),min(log10(T)),max(log10(T))]);
title(['$$\log_{10}\left(\Vert\hat{W_',num2str(r),'}\hat{W_',num2str(r),'}^T - \hat{U_',num2str(r),'}\hat{U_',num2str(r),'}^T\Vert_2\right)$$'],'Interpreter','latex');
end