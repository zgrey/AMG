% Model Adjoint Sensitivities of SCALED QiQi airfoil database
% Zach Grey - 08/26/2019
clc; close all; clearvars; rng(42);

addpath ~/AMG/Euclidean_tools/Euclidean_Shapes/;
cd ~/RESEARCH/SHDP/;

plt_flg = 1; msh_flg = 1; plt_modes = 10;
% things to modify
Nldmk = 1000; Nrnd = 1000; rPGA = 1000; tub = 1; tlb = -1;

%% read nominal points from file
dirstruct  = dir('~/airfoil_coordinate_database/coordinates/'); Naf = size(dirstruct,1)-2; 
% datapath = '/media/zgrey/46AFC5285FA7ADF9/AMG_DATA/Cd_adj'; dir_adj  = dir([datapath,'/adjoints']); Naf = size(dir_adj,1)-2;

%% transform airfoils using consistent affine map
P = zeros(Nldmk,2,Naf); b = zeros(2,1,Naf); i_fail = zeros(Naf,1); U = zeros(2*Nldmk,Naf);
M = zeros(2,2,Naf); Minv = zeros(2,2,Naf); n = zeros(Naf,1); not_tan = zeros(Naf,1);
for i=1:Naf
    clc; disp([num2str(i),'/',num2str(Naf),' Emb. Rep. shapes Complete']);
    % obtain nominal landmarks
    try 
        P0 = importairfoil([dirstruct(2+i).folder,'/',dirstruct(2+i).name]);
    catch
        disp(['Could not read: ',dirstruct(2+i).name])
        i_fail(i) = 1;
    end

    % use previously sucessfull meshes
%     data = csvread([dir_adj(2+i).folder,'/',dir_adj(2+i).name],1,0);
%     data_new = [data(1,:); data(3:end,:)];
%     P0 = [data_new(:,10) , data_new(:,11)];  

    % number of points
    n(i) = size(P0,1);

    %% compute embedding representation
    emb = embrep3(P0,Nldmk,'uni');
    % transform resampled points
    [Pemb,M(:,:,i),b(:,:,i),Minv(:,:,i)] = affine_trans(emb.pts,'LA'); P(:,:,i) = 1/sqrt(Nldmk-1)*Pemb;
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

% build section to original average scales (using average... because why not...)
proj_inv = sqrt(Nldmk - 1)*mean(Minv,3); 

% build section to original scales using Karcher mean
% disp('Computing Karcher mean for local section...');
% V = proj_inv; iter = 0; muM = eye(2);
% while norm(V,'fro') >= 1e-8 && iter <= 5
%     tic; sum_log = 0;
%     for i=1:Naf
%         sum_log = logm(M(:,:,i)*V) + sum_log;
%     end
%     V = V - 1/Naf*V*sum_log;
%     muM = Gr_exp(1,muM,V);
%     disp(['||V||_f = ',num2str(norm(V,'fro')),' ... ',num2str(toc),' sec.']);
%     iter = iter + 1;
% end
% muP0 = muP / muM';

%% PGA
disp('Computing PGA...')
for i=1:Naf
    Ui = Gr_log(muP,P(:,:,i));
    U(:,i) = reshape(Ui,2*Nldmk,1);
end
[PGA,PGAeigs] = svd(1/sqrt(Naf)*U);
% Uniform
trnd = tlb + (tub-tlb).*rand(Nrnd,rPGA);
% Gaussian
% trnd = randn(Nrnd,rPGA);

%% Plot PGA modes
fig = figure; PGAW = zeros(Nldmk,2,rPGA); PGA_fwd = zeros(Nldmk,2,rPGA);
sub = 10;
for i=1:plt_modes
    PGAW(:,:,i) = reshape((PGA(:,i)*PGAeigs(i,i)),Nldmk,2);
    PGA_fwd(:,:,i) = Gr_exp(1,muP,PGAW(:,:,i));
    subplot(2,5,i), bndplot(fig,muP(:,1),muP(:,2),sqrt(PGAW(:,1,i).^2 + PGAW(:,2,i).^2)); axis equal; hold on;
    subplot(2,5,i), plot(muP(:,1),muP(:,2),'w','linewidth',2); 
    subplot(2,5,i), plot(PGA_fwd(:,1,i),PGA_fwd(:,2,i),'k','linewidth',1);
    subplot(2,5,i), quiver(muP(1:sub:end,1),muP(1:sub:end,2),PGAW(1:sub:end,1,i),PGAW(1:sub:end,2,i),1,'k');
    fig.CurrentAxes.Visible = 'off';
end

%% plot and mesh airfoils
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

        % print figure
        print(fig,'-dpng',['./figs_PGA_airfoils/airfoil_',num2str(i),'.png'])
        % clear plot
        delete(h)
        
    end
    % write coordinates to *.dat file
    dlmwrite('airfoil.dat',[Pts(:,1,i) + 0.5,Pts(:,2,i), zeros(Nldmk,1)],'delimiter','\t');
    
    if msh_flg == 1
        % Mesh airfoil using gmsh 2.10.1
        try 
        fprintf('calling dat2gmsh\n');
        !python dat2gmsh.py airfoil.dat
        %     !python dat2gmsh_v2.py airfoil.dat
        fprintf('calling gmsh\n');
        % may require sudo for encrypted hard drives
        !gmsh airfoil.dat.geo -2 -o airfoil.mesh
        % refine mesh
%         !gmsh gmsh_refine -
        % save Gmsh geometry and mesh for visualization
        copyfile('./airfoil.mesh',['./PGA_gmsh_meshes/airfoil_',num2str(i),'.mesh']);
        copyfile('./airfoil.dat.geo',['./PGA_gmsh_geos/airfoil_',num2str(i),'.geo']);

        %% Convert Mesh
        tic;
        fprintf('mesh %i of %i pre-processing/conversion to .su2 format...\n',i,Nrnd);
        % Run DARPA EQUiPS SEQUOIA team codes:
        meshGMF = ReadGMF('airfoil.mesh');
        meshGMF = meshPrepro(meshGMF);
        meshSU2 = convertGMFtoSU2(meshGMF);
        WriteSU2_airfoil(meshSU2, ['./meshes/airfoil_',num2str(i),'.su2']);

        tocp = toc;
        % print time stats
        disp('SUCCESS!')
        disp(['Finished in... ',num2str(toc),' sec | Remaining time ~',num2str((Nrnd-i)*(toc+tocp)/120),' min']);
        catch
            disp('FAILED!')
            mesh_fail(i) = 1;
        end
    end
end
close all;
save('~/RESEARCH/AMG_DATA/qiqi_PGA_meshes.mat');