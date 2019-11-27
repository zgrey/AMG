% Grassmannian AMG post-processing
% Zach Grey - 08/26/2019
clc; close all; clearvars;

addpath ~/AMG/Euclidean_tools/Euclidean_Shapes/
addpath ~/RESEARCH/SHDP/

% datapath = '~/RESEARCH/SHDP';
% datapath = '/media/zgrey/46AFC5285FA7ADF9/AMG_DATA/Cd_adj'; QOI  = 3;
% datapath = '~/RESEARCH/AMG_DATA/Cd_adj'; QOI = 3;
datapath = '~/RESEARCH/AMG_DATA/PGA_samples/Cd_adj'; QOI = 3;
% datapath = '/media/zgrey/46AFC5285FA7ADF9/AMG_DATA/Cl_adj'; QOI = 2;

load([datapath,'/qiqi_PGA_meshes.mat'],'Minv_avg','PGA','rPGA');

% sweepdata = '/media/zgrey/46AFC5285FA7ADF9/AMG_DATA/Cd_adj/karcher';
% sweepdata = '~/RESEARCH/AMG_DATA/Cd_adj/karcher';
% sweepdata = '~/RESEARCH/AMG_DATA/PGA_samples/Cd_adj/karcher';

sweepdata = '~/RESEARCH/AMG_DATA/PGA_samples/Cd_adj/karcher/PGA_basis';

% number of resampled points
N = 1000; % 1000 achieved BEST results (visually)
% plots
plt_flg = 0; set(0,'defaulttextInterpreter','latex')
% meshes
msh_flg = 0; 
% differencing step size
h = 1e-4; % 1e-4 achieved BEST results (visually)
% smoothing
smth = 1;

% pick ranks
% rAS = 1000; rPGA = 4;
rAS = 400;

% pick central design
cnt_dsn = 'karcher';
% cnt_dsn = 'drag_opt_0AOA';
% cnt_dsn = 'drag_opt_1.25AOA';
% cnt_dsn = 'drag_opt_-1.25AOA';
% cnt_dsn = 'lift_opt';

%% adjoint smoothing bump function
if smth == 1
    s = 100; c = 0.05; bump = @(t) exp(s)./( (exp(c*s) + exp(s*t)).*(exp(c*s) + exp(s*(1-t))) );
elseif smth == 2
    
else
    bump = @(t) 1;
end

%% read adjoint information from file
dir_adj  = dir([datapath,'/adjoints']);
dir_force  = dir([datapath,'/iforces']);
Naf = size(dir_adj,1)-2;
H = zeros(N,2,Naf); Hproj = zeros(N,2,Naf); P = zeros(N,2,Naf);
b = zeros(2,1,Naf); M = zeros(2,2,Naf); Minv = zeros(2,2,Naf); n = zeros(Naf,1); 
not_tan = zeros(Naf,1); tan_err = zeros(Naf,1); Iadj = zeros(Naf,1);
for i=1:Naf
    clc; disp([num2str(i),'/',num2str(Naf),' Shapes Complete']);
    % read adjoints
    data = csvread([dir_adj(2+i).folder,'/',dir_adj(2+i).name],1,0);
    Idum = strfind(dir_adj(2+i).name,'_'); 
    
    % first naming convention
%     Iadj(i) = str2double(dir_adj(2+i).name(Idum(1)+1:Idum(2)-1));
    % simple naming convention
    Idum2 = strfind(dir_adj(2+i).name,'.'); Iadj(i) = str2double(dir_adj(2+i).name(Idum(1)+1:Idum2-1));
    % remove known LE and TE points reported as first two lines
    data_new = data(3:end,:);  
    % surface sensitivities
    surf_sens = data_new(:,9);
    % corresponding landmarks
    P0 = [data_new(:,10) , data_new(:,11)];  
    % number of points
    n(i) = size(P0,1);
     
    %% compute embedding representation
    emb = embrep3(P0,N,'uni');
    % transform resampled points
    [Pemb,M(:,:,i),b(:,:,i),Minv(:,:,i)] = affine_trans(emb.pts,'LA'); P(:,:,i) = 1/sqrt(N-1)*Pemb;
            
    %% compute perturbed embedding representation
    % sensitivity using x-sens and y-sens
    % THESE PERTURBATIONS ARE NOT NORMAL TO THE SURFACE?????? (SU2 bug?)
%     P_sens = bump(emb.nom.t).*[data_new(:,6) , data_new(:,7)];
    % compute pertrubations using surface sensitivity (SU2 assumes inward normals)
    P_sens = -bump(emb.nom.t).*surf_sens.*emb.nom.nml;

    % perturb physical scale landmarks
    % forward and backward perturbations 
    P_pert_fwd = P0 + h*P_sens; P_pert_bwd = P0 - h*P_sens;
    opts.AffineTrans = 'LA';
    emb_pert_bwd = embrep3(P_pert_bwd,N,'uni',opts);
    emb_pert_fwd = embrep3(P_pert_fwd,N,'uni',opts);
    % Affine-standardize perturbations (BEST!!!!!! same as simple SVD until the 5th airfoil????)
%     P_h_fwd = affine_trans(emb_pert_fwd.pts,'LA'); P_h_fwd = 1/sqrt(N-1)*P_h_fwd;
%     P_h_bwd = affine_trans(emb_pert_bwd.pts,'LA'); P_h_bwd = 1/sqrt(N-1)*P_h_bwd;
    % Affine-standardize perturbations using simple SVD
%     [P_h_fwd,~] = svd(emb_pert_fwd.pts - repmat(mean(emb_pert_fwd.pts,1),N,1),'econ');
%     [P_h_bwd,~] = svd(emb_pert_bwd.pts - repmat(mean(emb_pert_bwd.pts,1),N,1),'econ');
    % Affine-standardize with unperturbed transform 
%     [P_h_fwd,~] = svd(emb_pert_fwd.pts*M(:,:,i)' + repmat(b(:,:,i)',N,1),'econ');
%     [P_h_bwd,~] = svd(emb_pert_fwd.pts*M(:,:,i)' + repmat(b(:,:,i)',N,1),'econ');
    % transform perturbed points using unperturbed transform [Absil et al., Prop. 3.1]
    P_h_bwd = 1/sqrt(N-1)*(emb_pert_bwd.pts*M(:,:,i)' + repmat(b(:,:,i)',N,1));
    P_h_fwd = 1/sqrt(N-1)*(emb_pert_fwd.pts*M(:,:,i)' + repmat(b(:,:,i)',N,1));        
    
    %% compute shape sensitivity as Grassmannian gradient  
    % Taylor series approximation of gradient
    % LOOKS GOOD, WORTH RUNNING
%     H(:,:,i) = 1/h*Gr_log(P(:,:,i),P_h_bwd);
    % BEST RESULT
    H(:,:,i) = 1/h*Gr_log(P(:,:,i),P_h_fwd);
    % STRANGER LOOKING AMG PERTURBATIONS
%     H(:,:,i) = 1/(2*h)*(Gr_log(P(:,:,i),P_h_fwd) - Gr_log(P(:,:,i),P_h_bwd));
    % projection based gradient (THESE APPEAR QUITE NOISY, particularly at the LE)
     P_pert_spl1 = csape(emb.nom.t,P_sens(:,1),'periodic'); P_pert_spl2 = csape(emb.nom.t,P_sens(:,2),'periodic');  
     Hproj(:,:,i) = (eye(N) - P(:,:,i)*P(:,:,i)')*[bump(emb.t).*ppval(P_pert_spl1,emb.t) bump(emb.t).*ppval(P_pert_spl2,emb.t)];
%     H(:,:,i) = Hproj(:,:,i);
     
    tan_err(i) = norm((eye(N) - P(:,:,i)*P(:,:,i)')*H(:,:,i) - H(:,:,i),'fro');
    if tan_err(i) > 1e-8
        disp('WARNING: Some tangent vectors are not elements of tangent space - see not_tan')
        not_tan(i) = 1;
    end
    %% Visualize airfoils
    if plt_flg == 1
        fig = figure;
        % plot surface sensitivity
        bndplot(fig,emb.nom.pts(:,1),emb.nom.pts(:,2),log(sqrt(sum(P_sens.^2,2))));
        h1 = plot(emb.pts(:,1),emb.pts(:,2),'w','linewidth',2); axis equal; hold on; 
%         subplot(2,1,1), h3 = scatter(P0(:,1),P0(:,2),20,'filled','cdata',sqrt(sum(P_sens.^2,2))); colorbar; 
        h2 = quiver(emb.nom.pts(:,1),emb.nom.pts(:,2),P_sens(:,1),P_sens(:,2),1,'k','linewidth',2);
%         scatter(emb.nom.pts(:,1),emb.nom.pts(:,2),'k.')
        fig.CurrentAxes.Visible = 'off';
%         subplot(2,1,1), h1p = plot(emb_pert_fwd.pts(:,1),emb_pert_fwd.pts(:,2),'k','linewidth',1); axis equal; hold on;
%         subplot(2,1,1), h3pbwd = scatter(P_pert_bwd(:,1),P_pert_bwd(:,2),15,'k.');
%         subplot(2,1,1), h3pfwd = scatter(P_pert_fwd(:,1),P_pert_fwd(:,2),25,'k.');
        
        % plot adjoint components
        subplot(1,2,2), plot(emb.nom.t,P_sens(:,1),'linewidth',2); hold on;
        subplot(1,2,2), plot(emb.nom.t,P_sens(:,2),'linewidth',2);
        
        fig = figure; hviz = h; sub = 1;
        % Taylor series gradient
        plot(P(:,1,i),P(:,2,i),'linewidth',2); axis equal; hold on;
        plot(P_h_fwd(:,1),P_h_fwd(:,2),'k');
        quiver(P(1:sub:end,1,i),P(1:sub:end,2,i),H(1:sub:end,1,i),H(1:sub:end,2,i),1,'k');
        fig.CurrentAxes.Visible = 'off';
        % projected gradient
        fig = figure;
        plot(P(:,1,i),P(:,2,i),'linewidth',2); axis equal; hold on;
        quiver(P(1:sub:end,1,i),P(1:sub:end,2,i),Hproj(1:sub:end,1,i),Hproj(1:sub:end,2,i),1,'k')
        plot(P(:,1,i) + hviz*Hproj(:,1,i), P(:,2,i) + hviz*Hproj(:,2,i),'k')
        fig.CurrentAxes.Visible = 'off';
%           h6 = annotation('textbox', [0 0.9 1 0.1], ...
%             'String', ['Directory index i=',num2str(i),', ',dir_adj(2+i).name], ...
%             'EdgeColor', 'none', ...
%             'HorizontalAlignment', 'center');
        close all;

    end
end

if sum(not_tan) ~= 0
    disp('WARNING: Some tangent vectors are not elements of tangent space - see not_tan')
end
disp(['Average projection error of gradients = ',num2str(mean(tan_err))])

%% Central tangent space construction
if strcmp(cnt_dsn,'karcher')
    %% compute Kareshape(PTH0(:,2,:),N,Naf)rcher mean
    muP = P(:,:,1); V = ones(N,2); Log_P = zeros(N,2,Naf); iter = 0;
    disp('Computing Karcher mean...')
    disp('-----||V||_f History-----')
    while norm(V,'fro') >= 1e-8 && iter <= 5
        tic;
        for i=1:Naf
            if iter == 0 && i == 1
                Log_P(:,:,1) = zeros(N,2);
            else
                Log_P(:,:,i) = Gr_log(muP,P(:,:,i));
            end
        end
        V = 1/Naf*sum(Log_P,3);
        muP = Gr_exp(1,muP,V);
        disp(['||V||_f = ',num2str(norm(V,'fro')),' ... ',num2str(toc),' sec.']);
        iter = iter + 1;
    end

    % convert to original average scales (using average... because why not...)
    if ~exist('Minv_avg','var')
        disp('Computing inverse affine transformation for local section...');
        Minv_avg = sqrt(N - 1)*mean(Minv,3);
    end
    muP0 = muP*Minv_avg'; 

elseif strcmp(cnt_dsn,'drag_opt_0AOA')
    opt = csvread('/media/zgrey/46AFC5285FA7ADF9/AMG_DATA/SU2Opt_Cd_0AOA_surface_data.csv',1,0);
    muP0 = opt(:,10:11);
    muemb = embrep3(muP0,N,'uni');
    muPemb = affine_trans(muemb.pts,'LA'); 
    muP = 1/sqrt(N-1)*muPemb;
    Minv_avg = sqrt(N-1)*muemb.TF.Minv;
elseif strcmp(cnt_dsn,'drag_opt_1.25AOA')
    opt = csvread('/media/zgrey/46AFC5285FA7ADF9/AMG_DATA/SU2Opt_Cd_1.25AOA_surface_data.csv',1,0);
    muP0 = opt(:,10:11);
    muemb = embrep3(muP0,N,'uni');
    muPemb = affine_trans(muemb.pts,'LA'); 
    muP = 1/sqrt(N-1)*muPemb;
    Minv_avg = sqrt(N-1)*muemb.TF.Minv;
elseif strcmp(cnt_dsn,'drag_opt_-1.25AOA')
    opt = csvread('/media/zgrey/46AFC5285FA7ADF9/AMG_DATA/SU2Opt_Cd_-1.25AOA_surface_data.csv',1,0);
    muP0 = opt(:,10:11);
    muemb = embrep3(muP0,N,'uni');
    muPemb = affine_trans(muemb.pts,'LA'); 
    muP = 1/sqrt(N-1)*muPemb;
    Minv_avg = sqrt(N-1)*muemb.TF.Minv;
end

fig = figure;
subplot(1,2,1), plot(muP0(:,1),muP0(:,2),'k','linewidth',2); axis equal; hold on;
fig.CurrentAxes.Visible = 'off';
subplot(1,2,2), plot(muP(:,1),muP(:,2),'k','linewidth',2); axis equal; hold on;
fig.CurrentAxes.Visible = 'off';

%% Subspace computations (Parallel Trans. OPG, Embedded OPG, and PGA)
PTH0 = zeros(N,2,Naf); avgOPG = zeros(2*N); embOPG = zeros(2*N); U = zeros(2*N,Naf);
disp('Computing parallel translation of gradients...')
tic
for i=1:Naf
    % parallel translation of gradients
    PTH0(:,:,i) = Gr_parallel_trans(1,P(:,:,i),Gr_log(P(:,:,i),muP),H(:,:,i));
    PTgrad   = [reshape(PTH0(:,1,i),N,1);...
              reshape(PTH0(:,2,i),N,1)];
    avgOPG = PTgrad*PTgrad' + avgOPG;
    
    % embedded outer-product
    embgrad = [reshape(H(:,1,i),N,1);...
              reshape(H(:,2,i),N,1)];
    embOPG = embgrad*embgrad' + embOPG;
    
    % PGA
    Ui = Gr_log(muP,P(:,:,i));
    U(:,i) = reshape(Ui,2*N,1);
end
disp([num2str(toc),' sec.'])

% correlated SVD of Parallel Trans. OPG
[AS,Eigs] = svd(1/sqrt(Naf)*[reshape(PTH0(:,1,:),N,Naf);...
                   reshape(PTH0(:,2,:),N,Naf)]);
W = reshape(AS(:,1),N,2);
W1 = reshape(AS(:,1),N,2); W2 = reshape(AS(:,2),N,2);

% project to PGA basis
% rAS = rPGA;
% [AS,Eigs] = svd(1/sqrt(Naf)*PGA(:,1:rPGA)'*[reshape(PTH0(:,1,:),N,Naf);...
%                    reshape(PTH0(:,2,:),N,Naf)]);
% AS = PGA(:,1:rPGA)*AS;
% W = reshape(AS(:,1),N,2);
% W1 = reshape(AS(:,1),N,2); W2 = reshape(AS(:,2),N,2);

% embedded OPG [Mukherjee]
[embAS,embEigs] = svd(1/sqrt(Naf)*[reshape(H(:,1,:),N,Naf);...
                   reshape(H(:,2,:),N,Naf)]);
embW1 = reshape(embAS(:,1),N,2); embW2 = reshape(embAS(:,2),N,2);
% project to central tangent space
embW1 = (eye(N) - muP*muP')*embW1;
disp(['Norm difference of emb. and Par. Trans. basis = ',num2str(norm(embW1 - W,'fro'))])

% PGA
[PGA,PGAeigs] = svd(1/sqrt(Naf)*U);
PGAW = zeros(N,2,rPGA);
for i=1:rPGA, PGAW(:,:,i) = reshape(PGA(:,i),N,2); end

disp('Checking to see if AS basis is in horizontal space...')
disp(['Projection error = ',num2str(norm((eye(N) - muP*muP')*W - W,'fro'))])

% perturb small amount forward along AMG
AMG_fwd1 = Gr_exp(0.2,muP,W1); AMG_fwd2 = Gr_exp(0.1,muP,W2);
AMG0_fwd1 = AMG_fwd1*Minv_avg'; AMG0_fwd2 = AMG_fwd2*Minv_avg';
% perturb small amount forward along emb. OPG
embAMG_fwd1 = Gr_exp(0.2,muP,embW1); embAMG_fwd2 = Gr_exp(0.1,muP,embW2); 
embAMG0_fwd1 = embAMG_fwd1*Minv_avg'; embAMG0_fwd2 = embAMG_fwd2*Minv_avg';
% perturb small amount along PGA
PGA_fwd = zeros(N,2,rPGA);
for i=1:size(PGAW,3), PGA_fwd(:,:,i) = Gr_exp(0.1,muP,PGAW(:,:,i)); end

%% Visualize subspace computations
% visualize eigenvalues
figure; 
% subplot(1,2,1), scatter(1:Naf,diag(Eigs(1:Naf,1:Naf))); hold on;
subplot(1,2,1), scatter(1:rPGA,diag(Eigs(1:rPGA,1:rPGA))); hold on;
subplot(1,2,1), scatter(1:Naf,diag(embEigs(1:Naf,1:Naf)),15,'filled');
set(gca, 'YScale', 'log'); grid on; legend('AMG','Emb. OPG');
xlabel('$$index$$'); title('$$AMG \,\, Eigenvalues$$');
subplot(1,2,2), scatter(1:min([2*N,Naf]),diag(PGAeigs)); hold on;
set(gca, 'YScale', 'log'); grid on;
xlabel('$$index$$'); title('$$PGA \,\, Eigenvalues$$');

figure;
scatter(1:Naf,diag(Eigs(1:Naf,1:Naf))); hold on;
scatter(1:Naf,diag(embEigs(1:Naf,1:Naf)),10,'filled');
set(gca, 'YScale', 'log'); grid on;

figure;
scatter(1:min([2*N,Naf]),diag(PGAeigs)); hold on;
set(gca, 'YScale', 'log'); grid on;

% visualize eigenvectors
% quiver subset parameter
% first eigenvector
sub = 10; 
fig = figure;
% forward AMG
% subplot(2,1,1), bndplot(fig,muP(:,1),muP(:,2),sqrt(W(:,1).^2 + W(:,2).^2)); axis equal; hold on;
subplot(2,1,1), plot(muP(:,1),muP(:,2),'linewidth',2); axis equal; hold on;
subplot(2,1,1), plot(AMG_fwd1(:,1),AMG_fwd1(:,2),'k','linewidth',2);
subplot(2,1,1), quiver(muP(1:sub:end,1),muP(1:sub:end,2),W(1:sub:end,1),W(1:sub:end,2),1,'k');
subplot(2,1,1), embh = plot(embAMG_fwd1(:,1),embAMG_fwd1(:,2),'r--','linewidth',2);
fig.CurrentAxes.Visible = 'off';

% forward AMG over physical shape
% subplot(2,1,2), bndplot(fig,muP0(:,1),muP0(:,2),sqrt(W(:,1).^2 + W(:,2).^2)); axis equal; hold on; colorbar off; 
subplot(2,1,2), plot(muP0(:,1),muP0(:,2),'linewidth',2); axis equal; hold on;
subplot(2,1,2), plot(AMG0_fwd1(:,1),AMG0_fwd1(:,2),'k','linewidth',2);
subplot(2,1,2), plot(embAMG0_fwd1(:,1),embAMG0_fwd1(:,2),'--','color',embh.Color,'linewidth',2);
fig.CurrentAxes.Visible = 'off';

annotation('textbox', [0 0.9 1 0.1], ...
            'String', '$$1^{st}\,\,Eigenvector$$', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...Nt = 100;
            'interpreter','latex', ...
            'fontsize',16);
fig = figure;
% second eigenvector
% forward AMG
% subplot(2,1,1), bndplot(fig,muP(:,1),muP(:,2),sqrt(W2(:,1).^2 + W2(:,2).^2)); axis equal; hold on;
subplot(2,1,1), plot(muP(:,1),muP(:,2),'linewidth',2); axis equal; hold on;
subplot(2,1,1), plot(AMG_fwd2(:,1),AMG_fwd2(:,2),'k','linewidth',2);
subplot(2,1,1), quiver(muP(1:sub:end,1),muP(1:sub:end,2),W2(1:sub:end,1),W2(1:sub:end,2),1,'k');
subplot(2,1,1), embh = plot(embAMG_fwd2(:,1),embAMG_fwd2(:,2),'r--','linewidth',2);
fig.CurrentAxes.Visible = 'off';

% forward AMG over physical shape
% subplot(2,1,2), bndplot(fig,muP0(:,1),muP0(:,2),sqrt(W2(:,1).^2 + W2(:,2).^2)); axis equal; hold on; colorbar off;
subplot(2,1,2), plot(muP0(:,1),muP0(:,2),'linewidth',2); axis equal; hold on;
subplot(2,1,2), plot(AMG0_fwd2(:,1),AMG0_fwd2(:,2),'k','linewidth',2);
subplot(2,1,2), plot(embAMG0_fwd2(:,1),embAMG0_fwd2(:,2),'--','color',embh.Color,'linewidth',2);
fig.CurrentAxes.Visible = 'off';

annotation('textbox', [0 0.9 1 0.1], ...
            'String', '$$2^{nd}\,\,Eigenvector$$', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...Nt = 100;
            'interpreter','latex', ...
            'fontsize',16);
        
% PGA modes
fig = figure;
for i=1:4
%     subplot(2,2,i), bndplot(fig,muP(:,1),muP(:,2),sqrt(PGAW(:,1,i).^2 + PGAW(:,2,i).^2)); axis equal; hold on;
    subplot(2,2,i), plot(muP(:,1),muP(:,2),'linewidth',2); axis equal; hold on;
    subplot(2,2,i), plot(PGA_fwd(:,1,i),PGA_fwd(:,2,i),'k','linewidth',1);
    subplot(2,2,i), quiver(muP(1:sub:end,1),muP(1:sub:end,2),PGAW(1:sub:end,1,i),PGAW(1:sub:end,2,i),1,'k');
    fig.CurrentAxes.Visible = 'off';
end
annotation('textbox', [0 0.9 1 0.1], ...
            'String', '$$PGA\,\,Modes$$', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...Nt = 100;
            'interpreter','latex', ...
            'fontsize',16);
        
% Activity scores
ASscr = reshape(AS(:,1:rAS).^2*diag(Eigs(1:rAS,1:rAS)),N,2);

fig = figure;
% Parallel trans OPG
subplot(2,1,1), bndplot(fig,muP(:,1),muP(:,2),sqrt(ASscr(:,1).^2 + ASscr(:,2).^2)); 
axis equal; hold on;
subplot(2,1,1), plot(muP(:,1),muP(:,2),'w','linewidth',2); 
subplot(2,1,1), quiver(muP(1:sub:end,1),muP(1:sub:end,2),ASscr(1:sub:end,1),ASscr(1:sub:end,2),1,'k');
fig.CurrentAxes.Visible = 'off';

% Physical shape activity score
subplot(2,1,2), bndplot(fig,muP0(:,1),muP0(:,2),sqrt(ASscr(:,1).^2 + ASscr(:,2).^2));
axis equal; hold on; colorbar off;
subplot(2,1,2), plot(muP0(:,1),muP0(:,2),'w','linewidth',2); 
fig.CurrentAxes.Visible = 'off';

annotation('textbox', [0 0.9 1 0.1], ...
            'String', ['$$AMG\,\,activity\,\, Scores,\,\, r = $$',num2str(rAS)], ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'interpreter','latex', ...
            'fontsize',16);

% "project" to AMG
tproj = zeros(Naf,2);
disp('Optimizing to compute projections...')
tic;
for i=1:Naf
    tic;
    tproj(Iadj(i),1) = fminbnd(@(t) dGr_np(Gr_exp(t,muP,W1),P(:,:,i)),-2,2);
    tproj(Iadj(i),2) = fminbnd(@(t) dGr_np(Gr_exp(t,muP,W2),P(:,:,i)),-2,2);
    % simultaneous minimization
%     tproj(Iadj(i),:) = fmincon(@(t) dGr_np(Gr_exp(1,muP,...
%         reshape([reshape(W1,2*N,1), reshape(W2,2*N,1)]*t',N,2)),P(:,:,i)),[0,0],...
%         [],[],[],[],[-5,-5],[5,5]);
    
end
disp(['... finished in ',num2str(toc),' sec.']);

% build correlated sweep from concentrated measure
[Utprj,~] = svd(tproj','econ');
t2 = linspace(min(tproj*Utprj(:,1)),max(tproj*Utprj(:,1)),100);
t2sweep = t2'*Utprj(:,1)';

%% 1st principal AMG sweep
Nt = 100; 
% t = linspace(-0.5,0.5,Nt);
t = linspace(min(tproj(:,1)),max(tproj(:,1)),Nt);
% PICK SUBSPACE:
% dominant Parallel trans. AMG
WAMG = reshape(AS(:,1),N,2);
% correlated 2d Parallel trans. AMG
% WAMG = W1.*Utprj(1,1) + W2.*Utprj(2,1); t = t2;
% second dominant Parallel trans. AMG
% WAMG = reshape(AS(:,2),N,2); t = linspace(-0.2,0.2,Nt);
% dominant embedded OPG projection
% WAMG = embW;
% WAMG = PGAW(:,:,1);

AMG = zeros(N,2,length(t)); AMG0 = AMG; mesh_fail = zeros(Nt,1);

set(0,'defaulttextInterpreter','latex')
fig = figure; gifname = './AMG.gif';
for i=1:length(t)
    AMG(:,:,i)  = Gr_exp(t(i),muP,WAMG);
    
    % rescale and center
    AMG0(:,:,i) = AMG(:,:,i)*Minv_avg'; AMG0(:,:,i) = AMG0(:,:,i) - repmat(mean(AMG0(:,:,i),1),N,1);
    
    subplot(2,1,1), h1 = plot(AMG0(:,1,i),AMG0(:,2,i),'k','linewidth',2); hold on; axis equal;
    subplot(2,1,2), h2 = plot(AMG(:,1,i),AMG(:,2,i),'k','linewidth',2); hold on; axis equal;
    title(['t = ',num2str(t(i))]);
 
    % build gif
    figure(fig); frame = getframe(fig); 
    [A,map] = rgb2ind(frame2im(frame),256);
    if i == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.1);
    end
    delete([h1 h2]); 

    % mesh AMG
    if msh_flg == 1
    % shift to center in mesh a bit
    AMG0(:,:,i) = [AMG0(:,1,i) - min(AMG0(:,1,i)), AMG0(:,2,i)];
    % write coordinates to *.dat file
    dlmwrite('airfoil.dat',[AMG0(:,:,i), zeros(N,1)],'delimiter','\t');

    % Mesh airfoil using gmsh 2.10.1
    try 
        fprintf('calling dat2gmsh\n');
        !python dat2gmsh.py airfoil.dat
        fprintf('calling gmsh\n');
        !gmsh airfoil.dat.geo -2 -o airfoil.mesh
        % refine mesh
        %     !gmsh gmsh_refine -
        % save Gmsh geometry and mesh for visualization
        copyfile('./airfoil.mesh',['./gmsh_meshes/airfoil_',num2str(i),'.mesh']);
        copyfile('./airfoil.dat.geo',['./gmsh_geos/airfoil_',num2str(i),'.geo']);

        %% Convert Mesh
        tic;
        fprintf('mesh %i of %i pre-processing/conversion to .su2 format...\n',i,Nt);
        % Run DARPA EQUiPS SEQUOIA team codes:
        meshGMF = ReadGMF('airfoil.mesh');
        meshGMF = meshPrepro(meshGMF);
        meshSU2 = convertGMFtoSU2(meshGMF);
        WriteSU2_airfoil(meshSU2, ['./meshes/airfoil_',num2str(i),'.su2']);

        tocp = toc;
        % print time stats
        disp(['Finished in... ',num2str(toc),' sec | Remaining time ~',num2str((Nt-i)*(toc+tocp)/120),' min']);
    catch
        disp('MESH FAILED!')
        mesh_fail(i) = 1;
    end
    end
end

%% AMG shadows
% shadow over PGA
% W1 = reshape(PGA(:,1),N,2); W2 = reshape(PGA(:,2),N,2);

% read random function evaluations
if QOI  == 3
    % DRAG
    [Frnd,~,Iforces] = readsu2_forces(dir_force,QOI);
    ylbl = '$$C_d$$';
elseif QOI == 2
    % LIFT
    [Frnd,~,Iforces] = readsu2_forces(dir_force,QOI);
    ylbl = '$$C_{\ell}$$';
end

% scatter plots along AMG
figure; set(0,'defaulttextInterpreter','latex')
subplot(1,2,1), scatter(tproj(Iforces,1),Frnd(Iforces),50,'filled','cdata',Frnd(Iforces)); colorbar; hold on;
xlabel('$$ t $$'); ylabel('$$ f $$');
% if applicable, read sweep data
if exist([sweepdata,'/forces/airfoil_1.su2.dat'],'file')
    Fsweep = readsu2_forces(dir([sweepdata,'/forces']),QOI);
    %% determine warnings
    dirstruct = dir([sweepdata,'/WARNINGS']);
    N_warn = size(dirstruct,1)-2; 
    if N_warn ~= 0
        WARN = ones(N_warn,1);
        for i=1:N_warn
            ind1 = strfind(dirstruct(i+2).name,'_');
            ind2 = strfind(dirstruct(i+2).name,'.');
            WARN(i) = str2double(dirstruct(i+2).name(ind1+1:ind2-1));
        end
    end
    subplot(1,2,1), scatter(t,Fsweep,50,'filled','cdata',Fsweep);
    subplot(1,2,1), plot(t,Fsweep,'k','linewidth',2);
    
    % plot warnings
    scatter(t(WARN),Fsweep(WARN),200,'rx','linewidth',2)
end
subplot(1,2,2), scatter(tproj(Iforces,1),tproj(Iforces,2),50,'filled','cdata',Frnd(Iforces)); colorbar; hold on;
subplot(1,2,2), scatter(t,zeros(Nt,1),50,'filled','cdata',Fsweep);
% plot warnings
subplot(1,2,2), scatter(t(WARN),zeros(length(t(WARN)),1),200,'rx','linewidth',2)
%% save data
% close all;
clearvars sweepdata
save([datapath,'/',cnt_dsn,'/AMG_postproc.mat']);