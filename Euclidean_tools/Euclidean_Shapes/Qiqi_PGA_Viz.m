% Visualize Qiqi Airfoils

% Model Adjoint Sensitivities of SCALED QiQi airfoil database
% Zach Grey - 08/26/2019
clc; close all; clearvars; rng(42);

% Linux
HOME = '~/airfoil_coordinate_database';
% Windows
%HOME = 'C:\Users\zgrey\Documents\MATLAB\my_airfoil_coords';

% things to modify
plt_flg = 1;
Nldmk = 1000; Nrnd = 1000; rPGA = 10; tub = 1.5; tlb = -1.5;

%% read nominal points from file
dirstruct  = dir([HOME,'/coordinates/']); Naf = size(dirstruct,1)-2; 

%% transform airfoils using consistent affine map
P = zeros(Nldmk,2,Naf); b = zeros(2,1,Naf); i_fail = zeros(Naf,1); U = zeros(2*Nldmk,Naf);
M = zeros(2,2,Naf); Minv = zeros(2,2,Naf); n = zeros(Naf,1); not_tan = zeros(Naf,1); CovP0 = zeros(2,2,Naf);
for i=1:Naf
    clc; disp([num2str(i),'/',num2str(Naf),' Shapes Complete']);
    % obtain nominal landmarks
    try 
        P0 = importairfoil([dirstruct(2+i).folder,'/',dirstruct(2+i).name]);
    catch
        disp(['Could not read: ',dirstruct(2+i).name])
        i_fail(i) = 1; pause(1);
    end

    % number of points
    n(i) = size(P0,1);

    %% compute embedding representation
    emb = embrep3(P0,Nldmk,'uni');
    % transform resampled points
    [Pemb,M(:,:,i),b(:,:,i),Minv(:,:,i)] = affine_trans(emb.pts,'LA'); 
    P(:,:,i) = 1/sqrt(Nldmk-1)*Pemb;
    % compute convariance of recovered airfoil at consistent Nldmk
    CovP0(:,:,i) = 1/(Nldmk-1)*(emb.pts - repmat(mean(emb.pts,1),Nldmk,1))'*(emb.pts - repmat(mean(emb.pts,1),Nldmk,1));
end

%% compute Karcher mean of affine-standardized shapes
muP = P(:,:,1); V = ones(Nldmk,2); Log_P = zeros(Nldmk,2,Naf); iter = 0;
disp('Computing Karcher mean...')
disp('-----||V||_f History-----')
while norm(V,'fro') >= 1e-8 && iter <= 5
    tic;
    for i=1:Naf
        if iter == 0 && i == 1
            Log_P(:,:,1) = zeros(Nldmk,2);
        else
            Log_P(:,:,i) = Gr_log(muP,P(:,:,i));
        end
    end
    V = 1/Naf*sum(Log_P,3);
    muP = Gr_exp(1,muP,V);
    disp(['||V||_f = ',num2str(norm(V,'fro')),' ... ',num2str(toc),' sec.']);
    iter = iter + 1;
end

% Build constant local section of fiber bundle using arithmetic average
proj_inv = sqrt(Nldmk - 1)*mean(Minv,3);
% Compute physical-scale airfoil at Karcher mean
muP0 = muP*proj_inv';

%% PGA
disp('Computing PGA...')
for i=1:Naf
    Ui = Gr_log(muP,P(:,:,i));
    U(:,i) = reshape(Ui,2*Nldmk,1);
end
% compute PGA basis in central tangent space
[PGA,PGAeigs] = svd(1/sqrt(Naf)*U);
% compute corresponding Qiqi airfoil coordinates
t0 = U'*PGA(:,1:rPGA);
% generate i.i.d. Gaussian samples for generative model
trnd = randn(Nrnd,rPGA); 
% reassign the first samples first coord. to visualize
trnd(1,1) = -1;

% plot decay of eigenvalues
stem(diag(PGAeigs(1:rPGA,1:rPGA)));
set(gca, 'YScale', 'log'); grid on;
disp(['Soft threshold of variation: ',...
       num2str(sum(diag(PGAeigs(1:rPGA,1:rPGA)))/sum(diag(PGAeigs)))]);

%% Plot PGA coordinates
figure;
% first four coordinates
h = plotmatrix(t0(:,1:4)); hold on;

figure;
% all rPGA coordinates 
plotmatrix(t0); hold on;

%% Plot PGA modes
fig = figure; PGAW = zeros(Nldmk,2,rPGA); PGA_fwd = zeros(Nldmk,2,rPGA);
sub = 10;
for i=1:sub
    PGAW(:,:,i) = reshape((PGA(:,i)*PGAeigs(i,i)),Nldmk,2);
    PGA_fwd(:,:,i) = Gr_exp(1,muP,PGAW(:,:,i));
    subplot(2,5,i), bndplot(fig,muP(:,1),muP(:,2),sqrt(PGAW(:,1,i).^2 + PGAW(:,2,i).^2)); axis equal; hold on;
    subplot(2,5,i), plot(muP(:,1),muP(:,2),'w','linewidth',2); 
    subplot(2,5,i), plot(PGA_fwd(:,1,i),PGA_fwd(:,2,i),'k','linewidth',1);
    subplot(2,5,i), quiver(muP(1:sub:end,1),muP(1:sub:end,2),PGAW(1:sub:end,1,i),PGAW(1:sub:end,2,i),1,'k');
    fig.CurrentAxes.Visible = 'off';
end

%% Plot airfoils
mesh_fail = zeros(Naf,1); if plt_flg == 1, fig = figure; end
Pts = zeros(Nldmk,2,Naf); Gr_Pts = Pts;
for i=1:Nrnd
    % transform to scaled airfoil
    Vrnd = reshape((PGA(:,1:rPGA)*PGAeigs(1:rPGA,1:rPGA))*trnd(i,:)',Nldmk,2);
    Gr_Pts(:,:,i) = Gr_exp(1,muP,Vrnd);
    pts = Gr_Pts(:,:,i)*proj_inv';
    
    % shift to center in mesh a bit
    Pts(:,:,i) = pts;
    
    %% Visualize airfoils
    if plt_flg == 1
        % plot shape
        h = plot(Pts(:,1,i),Pts(:,2,i),'k','linewidth',2); axis equal; hold on;
        pause;
        delete(h);
    end
end