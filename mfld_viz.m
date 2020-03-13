clc; close all; clearvars;
%% Things to modify
% number of "data" points from manifold (for training)
N = 100;
% noise parameter for data
e = 0;
% number of intersecting quadratics
nq = 50;
% number of points per dimension for reevaluation
Nmfld = 100;
% size of ball for computing level-set intersection
int_ball = 0.05;
% multiple of furthest neighbor to size nearest neighbor for manifold-valued data
pct = 0.05;

%% model data-manifolds
% 1-manifold representations
p = 5; C = rand(p,3); t = linspace(-0.5,0.5,N);
% P0 = [ [cos(pi*t'), t'.^(1:p-1)]*C(:,1),...
%       [sin(pi*t').^2, t'.^(1:p-1)]*C(:,2),...
%       [cos(pi/2*t'), t'.^(1:p-1)]*C(:,3)]; 
P0 = [ cos(2*pi*t'), sin(2*pi*t'), 10*t'.^3];
% center
P0 = P0 - mean(P0,1);
% add noise (if any)
P = P0 + e*randn(N,1);
% visualize
fig = figure;
scatter3(P(:,1),P(:,2),P(:,3),50,'o','filled','r');
hold on; axis equal; grid on; alpha(0.5);
plot3(P0(:,1),P0(:,2),P0(:,3),'r--','linewidth',2)


% 2-manifold representations
% N = 100; S = [pi*rand(N,1), pi*rand(N,1)]; Nmfld = 100;
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

% quadratic(s)
cq = zeros(nchoosek(5,2),nq); G = cell(nq,1); H = cell(nq,1);
for i = 1:nq
    f = @(c,P) sum((P*[c(1) c(2) c(3);0 c(4) c(5);0 0 c(6)]').*P,2) + P*[c(7); c(8); c(9)] + c(10);

    % "training" or "learning"
    options = optimoptions('fminunc','Display','none');
    [cq(:,i),~] = fminunc(@(c) sum(f(c,P).^2),rand(nchoosek(5,2),1),options);
    % maps describing intrinsic properties of the trained zero-set
    A = [cq(1,i) cq(2,i) cq(3,i);0 cq(4,i) cq(5,i);0 0 cq(6,i)];
    a = [cq(7,i); cq(8,i); cq(9,i)];
    % gradient
    G{i} = @(P) P*(A' + A) + repmat(a',size(P,1),1); 
    % Hessian
    H{i} = A + A';
end

% sum of quadratic level-sets
[c3,~] = fminunc(@(c) sum(f(c(1:end/2),P).^2 + f(c(end/2+1:end),P).^2),10*rand(2*nchoosek(5,2),1),options);

% Fully connected neural net with ELU activation
% params
N3 = 1; N2 = 1;
strc = [3*ones(1,N3) 2*ones(1,N2) 1]; alph = 1;
% autoencoder
strc = [strc strc(end-1:-1:1)];
obj = @(c) sum(sum((P - NN_ELU(P, c(1:sum(strc.*[strc(2:end), 0])),...
                             c(sum(strc.*[strc(2:end), 0]) + 1:end),...
                             strc,alph) ).^2,2));
% prblm = createOptimProblem('fmincon','objective', obj,...
%                            'x0', 2*rand(sum(strc.*[strc(2:end), 0]) + sum(strc(2:end)),1)-1,...
%                            'lb',-100*ones(sum(strc.*[strc(2:end), 0]) + sum(strc(2:end)),1),...
%                            'ub',100*ones(sum(strc.*[strc(2:end), 0]) + sum(strc(2:end)),1));
% GS = GlobalSearch;
% [copt,LossOpt] = run(GS,prblm);
% vecW = copt(1:sum(strc.*[strc(2:end), 0]));
% vecb = copt(sum(strc.*[strc(2:end), 0]) + 1:end);

% random coefficients
vecW = randn(sum(strc.*[strc(2:end), 0]),1);
vecb = randn(sum(strc(2:end)),1);

%% resample domain
% high dimensional domain samples for patch surface
scl = 3;
[P1,P2,P3] = meshgrid(linspace(scl*min(P(:,1)),scl*max(P(:,1)),Nmfld),...
                      linspace(scl*min(P(:,2)),scl*max(P(:,2)),Nmfld),...
                      linspace(scl*min(P(:,3)),scl*max(P(:,3)),Nmfld));
PP = [reshape(P1,Nmfld^3,1), reshape(P2,Nmfld^3,1), reshape(P3,Nmfld^3,1)];

%% Visualize manifolds
% level-set of submersion 
FVnet = isosurface(P1, P2, P3,...
                   reshape(sum(NN_ELU(PP,vecW,vecb,...
                   strc(1:floor(end/2) +1),alph),2),...
                   Nmfld,Nmfld,Nmfld),0);
               
FVsum = isosurface(P1, P2, P3,reshape(f(c3(1:end/2),PP) + f(c3(end/2+1:end),PP),Nmfld,Nmfld,Nmfld),0);

FV = cell(nq,1);
for i=1:nq
    
    FV{i} = isosurface(P1, P2, P3,reshape(f(cq(:,i),PP),Nmfld,Nmfld,Nmfld),0);
    if ~isempty(FV{i}.vertices)
        % curvature of level set
        K = abs( sum((G{i}(FV{i}.vertices)*H{i}').*G{i}(FV{i}.vertices),2) - sum(G{i}(FV{i}.vertices).^2,2)*trace(H{i}) )...
            ./ (2*sum(G{i}(FV{i}.vertices).^2,2).^(3/2));

        % plot level-set
%         Ptch = patch(FV{i},'FaceVertexCData',log(K),'FaceColor','interp'); 
%         camlight; lighting phong; colorbar;
%         Ptch.FaceAlpha = 0.8; Ptch.EdgeColor = 'none';

        if i > 1
            D = pdist2(FVINT,FV{i}.vertices);
            FV_ind = max(D <= int_ball,[],1)'; FVINT= FV{i}.vertices(FV_ind,:);
        elseif i == 1
            FVINT = FV{i}.vertices;
        end
    end
end
                        
% intersection of quadratics
scatter3(FVINT(:,1),FVINT(:,2),FVINT(:,3),'k.');
% Ptch = patch(FVsum); 
% camlight; lighting phong; colorbar;
% Ptch.FaceAlpha = 0.5; Ptch.EdgeColor = 'none';

% nearest-neighbor over learned manifold
mfldD = pdist2(FVINT,P);
mfld_ind = max(mfldD <= pct*max(max(mfldD)),[],2)'; mfldPP = FVINT(mfld_ind,:);
scatter3(mfldPP(:,1),mfldPP(:,2),mfldPP(:,3),50,'ko','linewidth',2);

% Net submersion level-set (domain partition)
Ptch = patch(FVnet); 
camlight; lighting phong;
Ptch.FaceAlpha = 0.8; Ptch.EdgeColor = 'none';
alpha(0.5); colorbar;

% THIS ISN'T WORKING WELL (Bug?)
% neural-net manifold parameterization assuming strc = [3 2 1 2 3]
% size domain based on "learned" local chart
substrc = strc(1:(end + 1)/2);
subW = vecW(1:sum(substrc.*[substrc(2:end), 0]));
subb = vecb(1:sum(substrc(2:end)));
lclEuc = NN_ELU(P,subW,subb,substrc,alph);
% evaluate parameterization over compact domain
imW = vecW(7:end);
imb = vecb(4:end);
imstrc = substrc(end:-1:1);
NN_mfld = NN_ELU(linspace(min(lclEuc),max(lclEuc),100)', imW,imb,imstrc,alph);
figure; plot3(NN_mfld(:,1),NN_mfld(:,2),NN_mfld(:,3),'linewidth',2)

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