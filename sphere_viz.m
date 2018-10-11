% sphere visualizations
close all; clc;
%% The sphere
% mesh the entire sphere
nmesh = 50; [XX,YY,ZZ] = sphere(nmesh);
% parametrization of manifold (immersion)
X = @(S) [cos(S(:,1)).*sin(S(:,2)),sin(S(:,1)).*sin(S(:,2)),cos(S(:,2))]; 
% Jacobian of sphere parametrization
J = @(u,v) [-sin(u).*sin(v), cos(u).*cos(v); cos(u).*sin(v), sin(u).*cos(v); 0, -sin(v)];
% exponential map
Exp = @(t,Vt,P) kron(P,cos(t)) + kron(Vt./sqrt(sum(Vt.^2,2)),sin(t));
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
% w = 0.1; aa = 2*rand(m,1)-1; aa = aa/norm(aa); [A,~] = svd(aa); Func = @(XY) sum(sin(pi*XY*aa) + cos(pi/2*XY*aa),2) + w*sum((sin(pi*XY))*A(:,2:end),2); Grad = @(XY) kron(sum(pi*cos(pi*XY*aa) - pi/2*sin(pi/2*XY*aa),2),sum(aa,2)') + w*sum(pi*cos(pi*XY)*A(:,2:end),2)/(m-1); 
% highly nonlinear non-ridge
% aa = rand(m,1)-1; aa = aa/norm(aa); [A,~] = svd(aa); Func = @(XY) sum(sin(pi*XY*aa) + cos(pi/2*XY*aa),2) + sum((sin(pi*XY))*A(:,2:end),2); Grad = @(XY) kron(sum(pi*cos(pi*XY*aa) - pi/2*sin(pi/2*XY*aa),2),sum(aa,2)') + sum(pi*cos(pi*XY)*A(:,2:end),2)/(m-1);

%% Pushforward
% discretize manifold coordinates
N = 100; [u,v] = meshgrid(linspace(0,2*pi,sqrt(N)),linspace(0,pi,sqrt(N)));
S = reshape([u,v],N,2); P = X(S);
% domain tangent vectors:
% normalize points in domain
V = S./sqrt(sum(S.^2,2));
% or sample random directions
% V = 2*rand(N,2) - 1; V = V./sqrt(sum(V.^2,2));
% or a random constant direction
% V = 2*rand(1,2) -1; V = repmat(V/norm(V),N,1);
% or a constant coordinate direction
% V = repmat([0,1],N,1);

% compute tangent vectors to manifold (pushforward)
Vt = zeros(N,3); for ii=1:N, Vt(ii,:) = V(ii,:)*J(S(ii,1),S(ii,2))'; end
Vt = Vt./sqrt(sum(Vt.^2,2));

% Pushforward visualization
fig1 = figure; 
subplot(1,2,1), scatter(S(:,1),S(:,2),'filled','k'); hold on; axis equal; grid on; axis([-1,2*pi+1,-1,pi+1]);
subplot(1,2,1), quiver(S(:,1),S(:,2),V(:,1),V(:,2),1,'k');
subplot(1,2,2), mesh(XX,YY,ZZ,ones(size(ZZ))); hold on; axis equal; colormap gray; alpha(0.75);
subplot(1,2,2), scatter3(P(:,1),P(:,2),P(:,3),'filled','k');
subplot(1,2,2), quiver3(P(:,1),P(:,2),P(:,3),Vt(:,1),Vt(:,2),Vt(:,3),1,'k');
fig1.CurrentAxes.Visible = 'off';

%% Exponential map
N = 1; S = [1.4*pi, pi/2]; P = X(S); rng(2);
Vt = (2*rand(1,2) - 1)*J(S(1,1),S(1,2))'; Vt = Vt/norm(Vt);
% compute exponentials along tangent vector
k = 50; t = 2*linspace(0,1,k)'; geo = Exp(t,Vt,P);

% exponential map
fig2 = figure;
mesh(XX,YY,ZZ,ones(size(ZZ))); hold on; axis equal; colormap gray; alpha(0.75);
plot3(geo(:,1),geo(:,2),geo(:,3),'k','linewidth',2);
scatter3(P(1),P(2),P(3),50,'filled','g');
scatter3(geo(end,1),geo(end,2),geo(end,3),50,'filled','r');
quiver3(P(:,1),P(:,2),P(:,3),2*Vt(:,1),2*Vt(:,2),2*Vt(:,3),1,'b','linewidth',2);
fig2.CurrentAxes.Visible = 'off';

% Parallel translation approximation