%% Wrapper for example shape gradient evaluations
clc; close all; clearvars; rng(42);
load turbine_airfoil.mat;
addpath ./Euclidean_tools/Euclidean_Shapes/

%% parameter values
% t = 0;     % nominal
% t = 0.005; % average-area preserving
t = 0.01;  % area expanding
N = 500; p = 3/2; M = 100; q = 15;
A = (linspace(1,0,q).^p)'.*(2*rand(q,M)-1); B = (linspace(1,0,q).^p)'.*(2*rand(q,M)-1);

%% nominal shape
emb = embrep3(PP,N,'uni');
% nominal shape resampled at N points
PP0 = emb.pts;
nu0 = emb.nml;

%% turbine normalized field data
NN = 100;
[X,Y] = meshgrid(linspace(-2,2,NN),linspace(-1.5,1,NN));
XXX = [reshape(X,NN^2,1),reshape(Y,NN^2,1)];
XX = [field(:,1)-.1435, field(:,2)-0.065]; SS = 100*eye(2); mu = XXX;
FF = fit(XX,field(:,3),'cubicinterp'); w = FF(mu); I = find(~isnan(w));
w = w(I); mu = mu(I,:);
phi = rbf_surf(XX,w,mu,SS); 
phi(:,3) = phi(:,3)/max(phi(:,3))*max(field(:,3));

%% Loop shape perturbations
PP_nuT = zeros(N,2,M); nuT = zeros(N,2,M); GM = zeros(N,M);
for i=1:M

    % randomly perturbed shape
    emb = embpert3(emb,t,A(:,i),B(:,i),'expand');
    PP_nu = emb.pert.pts;
    nu = emb.pert.nml;
    kappa = emb.pert.curv;
    
    % pick curvature filter
    % radial basis filter
%     kappaRBF = rbf_surf(sarc,[kappa(1:1:end);kappa(1)],[sarc(1:1:end);1],50*N);
%     kappa = (kappaRBF(:,2)-min(kappaRBF(:,2)))/(max(kappaRBF(:,2))-min(kappaRBF(:,2)))...
%             *(max(kappa)-min(kappa))+min(kappa);
    % penalized
    kappa = sqrt(kappa.^2)*0.1;

    %% Generate rbf scalar field
    [phiS,~,GphiS] = rbf_surf(PP_nu,w,mu,SS);
    % DOUBLE CHECK THAT THIS QUANTITY CORRESPONDS TO LIFT
%     phiS(:,3) = (phiS(:,3)/max(phiS(:,3))*max(field(:,3)) ).*emb.pert.nml(:,2);
%     GphiS = GphiS.*[emb.pert.nml(:,2),emb.pert.nml(:,2)];
    
    % just pressure
    phiS(:,3) = (phiS(:,3)/max(phiS(:,3))*max(field(:,3)) );
    %% Gradient flow
    g  = sum(GphiS.*nu,2) + kappa.*phiS(:,3); 
    
    %% Save data
    PP_nuT(:,:,i) = PP_nu;
    nuT(:,:,i)    = nu;
    GM(:,i)       = g;
    
    if i==1
        %% Visualize shape characteristics
        skp = 1; % quiver can get congested, adjust to skip some entries
        fig = figure; hold on; axis equal;
        bndplot(fig,PP_nu(:,1),PP_nu(:,2),log(kappa));
        quiver(PP_nu(1:skp:end-1,1),PP_nu(1:skp:end-1,2),nu(1:skp:end-1,1),nu(1:skp:end-1,2));     
        plot(PP_nu(:,1),PP_nu(:,2),'k--');
        plot([PP(:,1);PP(1,1)],[PP(:,2);PP(1,2)],'k'); 
        ax = gca; ax.Visible = 'off';

        %% Visualize scalar field
        figure; hold on; axis equal; colorbar; title 'Scalar Field';
        scatter(XX(:,1),XX(:,2),'filled','cdata',phi(:,3));
        [phiS0,~,GphiS0] = rbf_surf(PP_nu,w,mu,SS);
        quiver(phiS0(1:skp:end,1),phiS0(1:skp:end,2),GphiS0(1:skp:end,1),GphiS0(1:skp:end,2))
        plot(PP_nu(:,1),PP_nu(:,2),'k','linewidth',2);
        ax = gca; ax.Visible = 'off';

        %% Visualize Gradient flow
        fig = figure; hold on; grid on; axis equal;
        bndplot(fig,PP_nu(:,1),PP_nu(:,2),g);
        ax = gca; ax.Visible = 'off';
    end

    clc; disp(['AGGREGATING DATA...',num2str(i/M*100),'% Complete']);
end

%% Subspace approximation
[U,Sig,D] = svd(GM); eps = -1;
%% Shape subspace plots
figure; scatter(1:length(diag(Sig)),diag(Sig).^2/Sig(1,1)^2,100,'filled'); grid on; set(gca,'yscale','log');
xlabel 'Index'
ylabel 'Eigenvalue'

% Perturb 4 shapes based on first singular vector
figure; hold on; axis equal
j = 1;
for i = [1,5,15,20]
    PP_nu = PP_nuT(:,:,i);
    nu = nuT(:,:,i);
    PP_eps = PP_nu + eps*nu.*U(:,1); %PP_eps = [PP_eps;PP_eps(1,:)];

    subplot(2,2,j), hold on; axis equal;
    subplot(2,2,j), plot([PP_nu(:,1);PP_nu(1,1)],[PP_nu(:,2);PP_nu(1,2)],'k'); 
    subplot(2,2,j), plot([PP_eps(:,1);PP_eps(1,1)],[PP_eps(:,2);PP_eps(1,2)],'k--');
    
    j=j+1;
end

% nominal perturbation
fig = figure; hold on; axis equal; PP_nu = PP0 + eps*nu0.*U(:,1);
bndplot(fig,PP0(:,1),PP0(:,2),-U(:,1));
plot(PP_nu(:,1),PP_nu(:,2),'k--');
plot([PP(:,1);PP(1,1)],[PP(:,2);PP(1,2)],'k');
ax = gca; ax.Visible = 'off';