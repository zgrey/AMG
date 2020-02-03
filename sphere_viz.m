% sphere visualizations
close all; clc;

%% The sphere
% mesh the entire sphere
nmesh = 50; [XX,YY,ZZ] = sphere(nmesh);
% parametrization of manifold (smooth immersion)
X = @(S) [cos(S(:,1)).*sin(S(:,2)),sin(S(:,1)).*sin(S(:,2)),cos(S(:,2))]; 
% Jacobian of sphere parametrization
J = @(u,v) [-sin(u).*sin(v), cos(u).*cos(v);...
             cos(u).*sin(v), sin(u).*cos(v); ...
             0, -sin(v)];
% exponential map
Exp = @(t,Vt,P) kron(P,cos(t)) + kron(Vt./sqrt(sum(Vt.^2,2)),sin(t));
% logarithmic map
Log = @(p,P) acos(P*p').*(P - P*p'.*repmat(p,size(P,1),1))./sqrt(sum((P - P*p'.*repmat(p,size(P,1),1)).^2,2));

%% Ambient map on the sphere
m = 3; XY = [reshape(XX,(nmesh+1)^2,1),reshape(YY,(nmesh+1)^2,1),reshape(ZZ,(nmesh+1)^2,1)]; rng(47);
% linear ambient function
a1 = 2*rand(m,1)-1; a1 = a1/norm(a1); Func = @(XY) XY*a1; Grad = @(XY) repmat(a1',size(XY,1),1);
% quadratic ambient ridge of rank(H) = r <= floor(m/2)
% r = 1; H = zeros(m); H(floor(m/2):floor(m/2)+r-1,floor(m/2):floor(m/2)+r-1) = eye(r); Func = @(X) sum((X*H).*X,2); Grad = @(X) X*H;
% highly nonlinear ridge
% a = 2*rand(m,1)-1; a = a/norm(a); Func = @(XY) sin(2*pi*XY*a) + cos(pi/2*XY*a); Grad = @(XY) kron(sum(pi*cos(pi*XY*a) - pi/2*sin(pi/2*XY*a),2),sum(a,2)');
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
speed = 1.75; k = 50; t = speed*linspace(0,1,k)'; geo = Exp(t,Vt,P);

% exponential map
fig2 = figure;
mesh(XX,YY,ZZ,ones(size(ZZ))); hold on; axis equal; colormap gray; alpha(1);
plot3(geo(:,1),geo(:,2),geo(:,3),'k','linewidth',2);
scatter3(P(1),P(2),P(3),80,'filled','g');
scatter3(geo(end,1),geo(end,2),geo(end,3),80,'filled','r');
quiver3(P(:,1),P(:,2),P(:,3),speed*Vt(:,1),speed*Vt(:,2),speed*Vt(:,3),1,'b','linewidth',2);
fig2.CurrentAxes.Visible = 'off';

%% Parallel translation approximation
N = 1; S = [1.4*pi, pi/2]; P0 = X(S); rng(2); k = 50; T = 1; t=T*linspace(0,1,k)';
Vt = (2*rand(1,2) - 1)*J(S(1,1),S(1,2))'; Vt = Vt/norm(Vt); P1 = Exp(T,Vt,P0); geo01 = Exp(t,Vt,P0);
Vp = (2*rand(1,2) - 1)*J(S(1,1),S(1,2))'; Vp = Vp/norm(Vp); v0 = Vp;

% central-differencing ladder
tau = 0.25; t=linspace(-1,1,50)';
P2    = Exp(tau,Vp,P0); P3 = Exp(-tau,Vp,P0); geo0 = Exp(tau*t,Vp,P0);
Vlog2 = Log(P1,P2); Vlog3 = Log(P1,P3);
diff = 1/(2*tau)*(Vlog2 - Vlog3);

% visualize ladder
fig3 = figure;
mesh(XX,YY,ZZ,ones(size(ZZ))); hold on; axis equal; colormap gray; alpha(1); fig3.CurrentAxes.Visible = 'off';
scatter3(P0(1),P0(2),P0(3),80,'filled','g');
scatter3(P1(1),P1(2),P1(3),80,'filled','r');
Pdiff = [P2;P3]; scatter3(Pdiff(:,1),Pdiff(:,2),Pdiff(:,3),80,'filled','b')
plot3(geo01(:,1),geo01(:,2),geo01(:,3),'k','linewidth',2);
plot3(geo0(:,1),geo0(:,2),geo0(:,3),'b','linewidth',2);
quiver3(P0(1),P0(2),P0(3),Vp(1),Vp(2),Vp(3),'g','linewidth',2);
quiver3([P1(1);P1(1)],[P1(2);P1(2)],[P1(3);P1(3)],[Vlog2(1);Vlog3(1)],[Vlog2(2);Vlog3(2)],[Vlog2(3);Vlog3(3)],'b','linewidth',2);
quiver3(P1(1),P1(2),P1(3),diff(1),diff(2),diff(3),'r','linewidth',2);

% visualize multiple rungs
fig4 = figure;
mesh(XX,YY,ZZ,ones(size(ZZ))); hold on; axis equal; colormap gray; alpha(1); fig4.CurrentAxes.Visible = 'off';
plot3(geo01(:,1),geo01(:,2),geo01(:,3),'k','linewidth',2);
nrung = 5; T = 1/nrung; P = Exp(linspace(0,1,nrung)',Vt,P0);
for i=1:nrung-1
    p0 = P(i,:); p1 = P(i+1,:);
    % central-differencing ladder
    tau = 0.1; t=linspace(-1,1,50)';
    P2    = Exp(tau,Vp,p0); P3 = Exp(-tau,Vp,p0); geo0 = Exp(tau*t,Vp,p0);
    Vlog2 = Log(p1,P2); Vlog3 = Log(p1,P3);
    diff = 1/(2*tau)*(Vlog2 - Vlog3);

    % visualize ladder
    scatter3(p0(1),p0(2),p0(3),80,'filled','g');
    scatter3(p1(1),p1(2),p1(3),80,'filled','r');
    Pdiff = [P2;P3]; scatter3(Pdiff(:,1),Pdiff(:,2),Pdiff(:,3),80,'filled','b')
    plot3(geo0(:,1),geo0(:,2),geo0(:,3),'b','linewidth',2);
    quiver3(p0(1),p0(2),p0(3),Vp(1),Vp(2),Vp(3),'g','linewidth',2);
    quiver3([p1(1);p1(1)],[p1(2);p1(2)],[p1(3);p1(3)],[Vlog2(1);Vlog3(1)],[Vlog2(2);Vlog3(2)],[Vlog2(3);Vlog3(3)],'b','linewidth',2);
    quiver3(p1(1),p1(2),p1(3),diff(1),diff(2),diff(3),'r','linewidth',2);

    % update vector transport
    Vp = diff;
end

% convergence of differencing ladder
N = 15; nrung=2.^linspace(0,N,N)'; err = ones(N,3); err_schilds = ones(N,1); err_ang = ones(N,1);
% compute angle to check for isometry
cost0 = v0*Vt'/(norm(v0)*norm(Vt));
for i=1:length(nrung)
    vp = diff_ladder(P0,P1,v0,Exp,Log,nrung(i));
    vp_s = schilds_ladder(P0,P1,v0,Exp,Log,nrung(i));
    err(i,:) = vp;
    err_schilds(i) = abs(-vp_s*Log(P1,P0)'/(norm(vp_s)*norm(Log(P1,P0))) - cost0);
    err_ang(i) = abs(-vp*Log(P1,P0)'/(norm(vp)*norm(Log(P1,P0))) - cost0);
end
err = err(1:end-1,:) - repmat(err(end,:),length(nrung)-1,1);
err = sqrt(err(:,1).^2 + err(:,2).^2 + err(:,3).^2);
figure; loglog(nrung,err_ang,'ko-','linewidth',2); grid on; hold on;
set(0,'defaulttextInterpreter','latex')
% loglog(nrung(1:end-1),err,'ko-');
% loglog(nrung,err_schilds,'bo-');
xlabel('number of rungs','fontsize',20); ylabel('error (preservation of angle)','fontsize',20);
% compute subspace distance convergence rate
M = [ones(length(nrung)-1,1) log10(nrung(1:end-1))]; cerr = M \ log10(err);
% coefficient of determination for convergence estimate
Rsq = 1 - sum((log10(err) - M*cerr).^2)/sum((log10(err) - mean(log10(err))).^2);
fprintf('Differencing ladder convergence rate 10^%f (R^2 = %f)\n',cerr(2),Rsq);

% Convergence of differencing schemes
h = linspace(eps,1,100);

%% Smooth 2-tensor study
N = 400; [u,v] = meshgrid(linspace(0,2*pi,sqrt(N)),linspace(0,pi,sqrt(N)));
S = reshape([u,v],N,2); P = X(S);
% random smooth curve
Ncurv = 200; 
t = linspace(0,2*pi,Ncurv)'; sc = [2*pi*(cos(t)+1), pi*(sin(t)+1)];
p0 = X(sc);

% compute random function and ambient gradient evaluations
F = Func(XY); Frnd = Func(P); Grnd = Grad(P);
% tangential gradient of ambient function
Gt = Grnd - bsxfun(@times,sum(Grnd.*P,2),P);

% compute a basis for each point tangent to the curve
% normalize points in domain
V = sc./sqrt(sum(sc.^2,2));
B1 = zeros(Ncurv,3); for ii=1:Ncurv, B1(ii,:) = V(ii,:)*J(sc(ii,1),sc(ii,2))'; end
B1 = B1./sqrt(sum(B1.^2,2)); B2 = cross(p0,B1,2);

% generate sphere viz
fig5 = figure;
% visualize map on sphere
surf(XX,YY,ZZ,reshape(F,(nmesh+1),(nmesh+1)),'FaceColor','interp','FaceLighting','gouraud','EdgeAlpha',0);
hold on; axis equal; alpha(0.25); fig5.CurrentAxes.Visible = 'off';
% scatter3(P(:,1),P(:,2),P(:,3),50,'filled','cdata',Frnd);
% quiver3(P(:,1),P(:,2),P(:,3),Gt(:,1),Gt(:,2),Gt(:,3),1,'k');
% visualize curve on sphere
plot3(p0(:,1),p0(:,2),p0(:,3),'k','linewidth',2);
% precondition AMG
U1 = zeros(Ncurv,3); U2 = zeros(Ncurv,3); err = zeros(Ncurv,1);
% precondition *.gif
filename = './smooth_2tensor.gif';
for i=1:Ncurv 
    % [central diff-ladder] 
    Vlog = zeros(N,3); Nrungs = 100;
    for ii=1:N
        Vlog(ii,:) = diff_ladder(P(ii,:),p0(i,:),Gt(ii,:),Exp,Log,Nrungs);
%         Vlog(ii,:) = schilds_ladder(P(ii,:),p0(i,:),Gt(ii,:),Exp,Log,Nrungs);
    end

    % SVD of parallel translation
    Ux = [B1(i,:);B2(i,:)];
    [U,D,~] = svd(1/sqrt(N)*Ux*(Vlog)',0); U = U'*Ux;
    U1(i,:) = U(1,:); U2(i,:) = U(2,:);
    
    % project ridge direction into each tangent space
    Proj_a = (eye(3) - p0(i,:)'*p0(i,:))*a1;
    % compute error
    err(i) = 1 - abs(U1(i,:)*Proj_a);
    % fix direction
    if U1(i,:)*Proj_a < 0, U1(i,:) = -U1(i,:); U2(i,:) = -U2(i,:); end
    
    Ptan = [p0(i,:) + B1(i,:); p0(i,:) + B2(i,:); p0(i,:) - B1(i,:); p0(i,:) - B2(i,:); p0(i,:) + B1(i,:);];
    % visualize tangent space at point
    h1 = scatter3(p0(i,1),p0(i,2),p0(i,3),50,'k','filled','MarkerEdgeColor','k','linewidth',2);
    h2 = plot3(Ptan(:,1),Ptan(:,2),Ptan(:,3),'k','linewidth',1);
    % plot AMG basis
    h3 = quiver3(p0(i,1),p0(i,2),p0(i,3),U1(i,1),U1(i,2),U1(i,3),1,'k','linewidth',2);
    h4 = quiver3(p0(i,1),p0(i,2),p0(i,3),U2(i,1),U2(i,2),U2(i,3),1,'k--','linewidth',2);
%     h4 = quiver3(p0(i,1)*ones(N,1),p0(i,2)*ones(N,1),p0(i,3)*ones(N,1),Vlog(:,1),Vlog(:,2),Vlog(:,3),'k--','linewidth',1);
    % plot projected ridge direction
    h5 = quiver3(p0(i,1),p0(i,2),p0(i,3),Proj_a(1),Proj_a(2),Proj_a(3),'r','linewidth',2);

    % rotate with tangent plane
%     view(p0(i,:));

    % build gif
    figure(fig5); frame = getframe(fig5);
    [A1,map] = rgb2ind(frame2im(frame),256);
    if i == 1
        imwrite(A1,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A1,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    
    delete([h1,h2,h3,h4,h5]);
end
