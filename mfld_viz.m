clc; close all; clearvars;
%% model manifolds
% 1-manifold representations
N = 100; e = 0; p = 5; Nmfld = 100;
C = rand(p,3); t = linspace(-0.5,0.5,N);
% P0 = [ [cos(pi*t'), t'.^(1:p-1)]*C(:,1),...
%       [sin(pi*t').^2, t'.^(1:p-1)]*C(:,2),...
%       [cos(pi/2*t'), t'.^(1:p-1)]*C(:,3)]; 
P0 = [ cos(2*pi*t'), sin(2*pi*t'), 10*t'.^3];
P = P0 + e*randn(N,1);
fig = figure;
scatter3(P(:,1),P(:,2),P(:,3),50,'o','filled','k');
hold on; axis equal; grid on; alpha(0.5);
plot3(P0(:,1),P0(:,2),P0(:,3),'r--','linewidth',2)

% 2-manifold representations
% N = 100; S = [2*pi*rand(N,1), pi/2*rand(N,1)]; Nmfld = 100;
% % sphere
% X = @(S) [cos(S(:,1)).*sin(S(:,2)),sin(S(:,1)).*sin(S(:,2)),cos(S(:,2))];
% % torus
% % X = @(S) 0.25*[(2 + cos(S(:,2))).*cos(S(:,1)), (2 + cos(S(:,2))).*sin(S(:,1)), sin(S(:,2))];
% P = X(S);
% fig = figure;
% [U,V] = meshgrid(linspace(0,2*pi,Nmfld),linspace(0,2*pi,Nmfld));
% M = X([reshape(U,Nmfld^2,1), reshape(V,Nmfld^2,1)]);
% XX = reshape(M(:,1),Nmfld,Nmfld); YY = reshape(M(:,2),Nmfld,Nmfld); ZZ = reshape(M(:,3),Nmfld,Nmfld);
% h = mesh(XX,YY,ZZ,ones(size(ZZ))); h.EdgeColor = [0.3,0.75,0.93];
% hold on; axis equal; alpha(0.5); 
% scatter3(P(:,1),P(:,2),P(:,3),50,'o','filled','markeredgecolor','k'); alpha(0.5);

%% manifold learning
% subspace
% [A1,D] = svd((P - mean(P,1))');
% [Y1,Y2] = meshgrid(linspace(-1,1,Nmfld),linspace(-1,1,Nmfld));
% Ptan = [mean(P,1) + A1(:,1)'; mean(P,1) + A1(:,2)'; mean(P,1) - A1(:,1)';...
%         mean(P,1) - A1(:,2)'; mean(P,1) + A1(:,1)'];
% plot3(Ptan(:,1),Ptan(:,2),Ptan(:,3),'k-','linewidth',2);
% quiver3(mean(P(:,1)),mean(P(:,2)),mean(P(:,3)),A(1,1),A(2,1),A(3,1),1,'k','linewidth',2)
% quiver3(mean(P(:,1)),mean(P(:,2)),mean(P(:,3)),A(1,2),A(2,2),A(3,2),1,'k','linewidth',2)
% alpha(0.5);

% graph of function
% V = @(P) [ones(size(P,1),1) P(:,1) P(:,2) P(:,1).*P(:,2) P(:,1).^2 P(:,2).^2];
% c = V(P) \ P(:,3);
% F = V([reshape(Y1,Nmfld^2,1), reshape(Y2,Nmfld^2,1)])*c;
% h = surf(Y1,Y2,reshape(F,Nmfld,Nmfld)); alpha(0.5);
% h.FaceAlpha = 0.5;
% h.EdgeColor = 'none';

% quadratic
f = @(c,P) sum((P*[c(1) c(2) c(3);0 c(4) c(5);0 0 c(6)]').*P,2) + P*[c(7); c(8); c(9)] + c(10);

% "training" or "learning"
[c1,fopt1] = fminunc(@(c) sum(f(c,P).^2),rand(nchoosek(5,2),1));
% maps describing intrinsic properties of the trained zero-set
A1 = [c1(1) c1(2) c1(3);0 c1(4) c1(5);0 0 c1(6)]; a1 = [c1(7); c1(8); c1(9)];
% gradient
G1 = @(P) P*(A1' + A1) + repmat(a1',size(P,1),1); 
% Hessian
H1 = A1 + A1';

% another quadratic
[c2,fopt2] = fminunc(@(c) sum(f(c,P).^2),10*rand(nchoosek(5,2),1));
A2 = [c2(1) c2(2) c2(3);0 c2(4) c2(5);0 0 c2(6)]; a2 = [c2(7); c2(8); c2(9)];
% gradient
G2 = @(P) P*(A2' + A2) + repmat(a2',size(P,1),1); 
% Hessian
H2 = A2 + A2';

% "intersection" of quadratics
[c3,~] = fminunc(@(c) sum(f(c(1:end/2),P).^2 + f(c(end/2+1:end),P).^2),10*rand(2*nchoosek(5,2),1));

% "encoder"
% params
d = 2; m = 3;
obj = @(c) sum(Enc_ELU(P,reshape([c(1:m^2*d - 2*m); zeros(2*m,1)],m,m,d),...
                         reshape([c(m^2*d - 2*m + 1:end); zeros(m-1,1)], m,d) ).^2);
prblm = createOptimProblem('fmincon','objective', obj,...
                        'x0', rand(m^2*d - 2*m + m*d - (m - 1),1));
GS = GlobalSearch;
[copt,~] = run(GS,prblm);

% build aggregate of "trained" affine functions
Wopt = reshape([copt(1:m^2*d - 2*m); zeros(2*m,1)],m,m,d);
bopt = reshape([copt(m^2*d - 2*m + 1:end); zeros(m-1,1)], m,d);

% "autoencoder" (Decoder o Encoder)
% shft = m^2*d - 2*m + m*d - (m - 1) +1;
% obj = @(c) sum(Dec_ELU(Enc_ELU(P,reshape([c(1:m^2*d - 2*m); zeros(2*m,1)],m,m,d),...
%                          reshape([c(m^2*d - 2*m + 1:end); zeros(m-1,1)], m,d)),...
%                          reshape([zeros(2*m,1); c(shft:shft + m^2*d - 2*m)],m,m,d),...
%                          reshape([zeros(m-1,1); c(shft + m^2*d - 2*m + 1:end)], m,d) ).^2);
% prblm = createOptimProblem('fmincon','objective', obj,...
%                         'x0', rand(2*(m^2*d - 2*m + m*d - (m - 1)),1));
% GS = GlobalSearch;
% [copt,~] = run(GS,prblm);

%% resample domain
% high dimensional domain samples for patch surface
scl = 3;
[P1,P2,P3] = meshgrid(scl*linspace(min(P(:,1)),max(P(:,1)),Nmfld),...
                      scl*linspace(min(P(:,2)),max(P(:,2)),Nmfld),...
                      scl*linspace(min(P(:,3)),max(P(:,3)),Nmfld));
% [P1,P2,P3] = meshgrid(linspace(-1,1,Nmfld),...
%                       linspace(-1,1,Nmfld),...
%                       linspace(-1,1,Nmfld));

%% Visualize shapes
PP = [reshape(P1,Nmfld^3,1), reshape(P2,Nmfld^3,1), reshape(P3,Nmfld^3,1)];
FV1 = isosurface(P1, P2, P3,reshape(f(c1,PP),Nmfld,Nmfld,Nmfld),0);
FV2 = isosurface(P1, P2, P3,reshape(f(c2,PP),Nmfld,Nmfld,Nmfld),0);
FV1cap2 = isosurface(P1, P2, P3,reshape(f(c3(1:end/2),PP) + f(c3(end/2+1:end),PP),Nmfld,Nmfld,Nmfld),0); 
FVnet = isosurface(P1, P2, P3,reshape(Enc_ELU(PP,Wopt,bopt),Nmfld,Nmfld,Nmfld),0);

% curvature of level set
K1 = abs( sum((G1(FV1.vertices)*H1').*G1(FV1.vertices),2) - sum(G1(FV1.vertices).^2,2)*trace(H1) )...
    ./ (2*sum(G1(FV1.vertices).^2,2).^(3/2));
K2 = abs( sum((G2(FV2.vertices)*H2').*G2(FV2.vertices),2) - sum(G2(FV2.vertices).^2,2)*trace(H2) )...
    ./ (2*sum(G2(FV2.vertices).^2,2).^(3/2));

% first quadratic
Ptch = patch(FV1,'FaceVertexCData',log(K1),'FaceColor','interp'); 
camlight; lighting phong; colorbar;
Ptch.FaceAlpha = 0.8; Ptch.EdgeColor = 'none';

% second quadratic
Ptch = patch(FV2,'FaceVertexCData',log(K2),'FaceColor','interp'); 
camlight; lighting phong; colorbar;
Ptch.FaceAlpha = 0.8; Ptch.EdgeColor = 'none';

% intersection of quadratics
% Ptch = patch(FV1cap2); 
% camlight; lighting phong; colorbar;
% Ptch.FaceAlpha = 0.5; Ptch.EdgeColor = 'none';

% Net
Ptch = patch(FVnet); 
camlight; lighting phong;
Ptch.FaceAlpha = 0.8; Ptch.EdgeColor = 'none';
% Dense visualization of Net response
% figure; scatter3(PP(1:70:end,1),PP(1:70:end,2),PP(1:70:end,3),25,'filled','cdata',Enc_ELU(PP(1:70:end,:),Wopt,bopt));
% alpha(0.5); colorbar;

%% Convergence visualization
% imgfig = figure;
% h = mesh(XX,YY,ZZ,ones(size(ZZ))); h.EdgeColor = [0.3,0.75,0.93];
% hold on; axis equal; alpha(0.5); imgfig.CurrentAxes.Visible = 'off';
% gifname = './img/mfld_learn.jpg';
% for N = 1
%     S = [2*pi*rand(N,1), pi*rand(N,1)];
%     P = X(S);
%     
%     % data
%     h = scatter3(P(:,1),P(:,2),P(:,3),50,'o','filled','markeredgecolor','k');
%     
%     [c,fopt] = fminunc(@(c) sum(f(c,P).^2),rand(nchoosek(5,2),1));
%     [P1,P2,P3] = meshgrid(2*linspace(-1,1,Nmfld),2*linspace(-1,1,Nmfld),2*linspace(-1,1,Nmfld));
%     PP = [reshape(P1,Nmfld^3,1), reshape(P2,Nmfld^3,1), reshape(P3,Nmfld^3,1)];
%     FV = isosurface(P1, P2, P3,reshape(f(c,PP),Nmfld,Nmfld,Nmfld),0);
%     Ptch = patch(FV);
%     Ptch.FaceAlpha = 0.5;
%     Ptch.EdgeColor = 'none';
%     Ptch.FaceColor = [0.3,0.75,0.93];
%     axis([-2,2,-2,2]);
%     
%     % build gif
%     figure(imgfig); frame = getframe(imgfig); 
%     [A,map] = rgb2ind(frame2im(frame),256);
%     imwrite(A,map,[gifname,'-',num2str(N),'.jpg'],'jpeg');
%         
%     delete([h,Ptch]);
%     
% end

%% Rotating axis movie
% get current figure
% fig = gcf;
% filename = 'mfld4.gif';
% theta = linspace(0,2*pi,100);
% [AZ,EL] = view;
% for i=1:length(theta)
%     view(AZ+theta(i)*180/pi,EL);
%     frame = getframe(fig);
%     [A1,map] = rgb2ind(frame2im(frame),256);
%     if i == 1
%         imwrite(A1,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);    
%     else
%         imwrite(A1,map,filename,'gif','WriteMode','append','DelayTime',0.1);
%     end
% end