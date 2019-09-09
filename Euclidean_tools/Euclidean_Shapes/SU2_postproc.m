% Grassmannian AMG post-processing
% Zach Grey - 08/26/2019
clc; close all; clearvars;

addpath ~/AMG/Euclidean_tools/Euclidean_Shapes/
addpath ~/RESEARCH/SHDP/

% datapath = '~/RESEARCH/SHDP/adjoints';
datapath = '/media/zgrey/46AFC5285FA7ADF9/AMG_DATA/disc_adj/adjoints';

% number of resampled points
N = 1000; plt_flg = 0; msh_flg = 1; h = 5e-4;

% smoothing bump function
r = 100; c = 0.05; bump = @(t) exp(r)./( (exp(c*r) + exp(r*t)).*(exp(c*r) + exp(r*(1-t))) );

%% read adjoint information from file
dirstruct  = dir(datapath);
Naf = size(dirstruct,1)-2;
H = zeros(N,2,Naf); Hproj = zeros(N,2,Naf); P = zeros(N,2,Naf);
b = zeros(2,1,Naf); Minv = zeros(2,2,Naf); n = zeros(Naf,1);
disp('Reading adjoint data...');
for i=1:Naf
    clc; disp([num2str(i),'/',num2str(Naf),' Shapes Complete']);
    data = csvread([dirstruct(2+i).folder,'/',dirstruct(2+i).name],1,0);
    % remove known LE and TE points reported as first two lines
    data_new = data(3:end,:);
%     data_new = [data(1,:); data(3:end,:)];
    
    surf_sens = data_new(:,9);
    P_sens = [data_new(:,6) , data_new(:,7)];
    P0 = [data_new(:,10) , data_new(:,11)];
    
    % number of points
    n(i) = size(P0,1);
    
    %% compute embedding representation
    emb = embrep3(P0,N,'uni');
    TF = emb.pts*emb.TF.M' + repmat(emb.TF.b',N,1);
    b(:,:,i) = emb.TF.b; Minv(:,:,i) = emb.TF.Minv;
    
    %% compute perturbed embedding representation
    % perturbation using x-sens and y-sens
%     P_pert = P0 + h*P_sens;
    % perturbation using surface sensitivity assuming inward normals
%     P_sens = [surf_sens.*emb.nom.nml(:,1) surf_sens.*emb.nom.nml(:,2)]; P_pert = P0 - h*P_sens;
    
    % perturbation using penalized x-sens and y-sens
%     P_sens = bump(emb.nom.t).*P_sens; P_pert = P0 + h*P_sens;
    % perturbation of penalized surface sensitivity
    P_sens = bump(emb.nom.t).*surf_sens.*emb.nom.nml; P_pert = P0 + h*P_sens;
    
    % perturb
    emb_pert = embrep3(P_pert,N,'uni');
    TF_sens = emb_pert.pts*emb_pert.TF.M' + repmat(emb_pert.TF.b',N,1);
    
    %% compute shape sensitivity as tangent vector
    % Taylor series approximation
    [P(:,:,i),~] = svd(TF,'econ'); [P_h,~] = svd(TF_sens,'econ');
    H(:,:,i) = 1/h*Gr_log(P(:,:,i),P_h);
    % projection based gradient
    P_pert_spl1 = csape(emb.nom.t,P_sens(:,1),'periodic'); P_pert_spl2 = csape(emb.nom.t,P_sens(:,2),'periodic');  
    Hproj(:,:,i) = (eye(N,N) - P(:,:,i)*P(:,:,i)')*[bump(emb.t).*ppval(P_pert_spl1,emb.t) bump(emb.t).*ppval(P_pert_spl2,emb.t)];
    
    %% Visualize airfoils
    if plt_flg == 1
        % plot continuous approximation of shape
        subplot(1,2,1), h1 = plot(emb.pts(:,1),emb.pts(:,2),'linewidth',2); axis equal; hold on;
%         subplot(1,2,1), h2 = quiver(emb.pts(:,1),emb.pts(:,2),emb.nml(:,1),emb.nml(:,2),'color',h1.Color);
        subplot(1,2,1), h1p = plot(emb_pert.pts(:,1),emb_pert.pts(:,2),'linewidth',2); axis equal; hold on;
        subplot(1,2,1), h2p = quiver(emb_pert.pts(:,1),emb_pert.pts(:,2),emb_pert.nml(:,1),emb_pert.nml(:,2),'color',h1p.Color);
        % plot original airfoil
        subplot(1,2,1), h3 = scatter(P0(:,1),P0(:,2),30,'filled','cdata',surf_sens);
        subplot(1,2,1), h3p = scatter(P_pert(:,1),P_pert(:,2),15,'k.');
        fig.CurrentAxes.Visible = 'off';
        % plot transformed coordinates
        subplot(1,2,2), h4 = plot(emb.TF.pts(:,1),emb.TF.pts(:,2),'linewidth',2); hold on;
        subplot(1,2,2), h7 = plot(TF(:,1),TF(:,2),'--','linewidth',2);
        subplot(1,2,2), h5 = scatter(emb.TF.pts(:,1),emb.TF.pts(:,2),15,'k');
        % plot transformed perturbed coordinates
        subplot(1,2,2), h7p = plot(TF_sens(:,1),TF_sens(:,2),'--','linewidth',2);
        subplot(1,2,2), h5p = scatter(emb_pert.TF.pts(:,1),emb_pert.TF.pts(:,2),15,'k.');
        
        figure; hviz = 0.01;
        % Taylor series gradient
        subplot(1,2,1), plot(P(:,1,i),P(:,2,i),'linewidth',2); axis equal; hold on;
        subplot(1,2,1), plot(P(:,1,i) + hviz*H(:,1,i), P(:,2,i) + hviz*H(:,2,i),'linewidth',2)
        % projected gradient
        subplot(1,2,2), plot(P(:,1,i),P(:,2,i),'linewidth',2); axis equal; hold on;
        subplot(1,2,2), plot(P(:,1,i) + hviz*Hproj(:,1,i), P(:,2,i) + hviz*Hproj(:,2,i),'linewidth',2)
        h6 = annotation('textbox', [0 0.9 1 0.1], ...
            'String', ['Directory index i=',num2str(i),', ',dirstruct(2+i).name], ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');[AS,Eig,inAS] = eigs(avgOPG);
        close all;

    end
end

%% compute Kareshape(PTH0(:,2,:),N,Naf)rcher mean
muP = P(:,:,1); V = ones(N,2); Log_P = zeros(N,2,Naf); iter = 0;
disp('Computing Karcher mean...')
disp('-----||V||_f History-----')
while norm(V,'fro') >= 1e-8 && iter < 1e4
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

% convert to original average scales (using avera0.25ge, because why not)
Minv_avg = sqrt(mean(n) - 1)*mean(Minv,3);
muP0 = muP*Minv_avg'; 

fig = figure;
subplot(1,2,1), 
subplot(1,2,1), plot(muP0(:,1),muP0(:,2),'linewidth',2); axis equal; hold on;
fig.CurrentAxes.Visible = 'off';
subplot(1,2,2), plot(muP(:,1),muP(:,2),'linewidth',2); axis equal; hold on;

%% Parallel translate gradients
PTH0 = zeros(N,2,Naf); avgOPG = zeros(N,N);
disp('Computing parallel translation of gradients...')
tic
for i=1:Naf
    PTH0(:,:,i) = Gr_parallel_trans(1,P(:,:,i),Gr_log(P(:,:,i),muP),H(:,:,i));
%     PTH0(:,:,i) = Gr_parallel_trans(1,P(:,:,i),Gr_log(P(:,:,i),muP),Hproj(:,:,i));
end
disp([num2str(toc),' sec.'])

% independent SVD
[AS1,Lambda1] = svd(reshape(PTH0(:,1,:),N,Naf));
[AS2,Lambda2] = svd(reshape(PTH0(:,2,:),N,Naf));
% aggregate important directions
W = [AS1(:,1) AS2(:,1)];

% correlated SVD
[AS,Lambda] = svd([reshape(PTH0(:,1,:),N,Naf);...
                   reshape(PTH0(:,2,:),N,Naf)]);
W = reshape(AS(:,1),N,2); W = W./repmat([norm(W(:,1)), norm(W(:,2))],N,1);

%% compute AMG
Nt = 100;
t = linspace(-0.3,0.3,Nt);
AMG = zeros(N,2,length(t)); AMG0 = AMG; mesh_fail = zeros(Nt,1);

set(0,'defaulttextInterpreter','latex')
fig = figure; gifname = './AMG.gif';
% plot(muP0(:,1),muP0(:,2),'--'); hold on;
% P0 = P(:,:,1)*Minv_avg'; plot(P0(:,1),P0(:,2),'-.');
for i=1:length(t)
    AMG(:,:,i)  = Gr_exp(t(i),muP,W);
    AMG0(:,:,i) = AMG(:,:,i)*Minv_avg';
    
    subplot(2,1,1), h1 = plot(AMG0(:,1,i),AMG0(:,2,i),'k','linewidth',2); hold on; axis equal;
    subplot(2,1,2), h2 = plot(AMG(:,1,i),AMG(:,2,i),'k','linewidth',2); hold on; axis equal;
    subplot(2,1,1), fig.CurrentAxes.Visible = 'off';
    subplot(2,1,2), fig.CurrentAxes.Visible = 'off';
    title(['t = ',num2str(t(i))]);
    
    % shape geodesic
%     [H,~,U,S,V] = Gr_log(P(:,:,1),muP);
%     Pgeo = Gr_exp(t(i),P(:,:,1),H)*Minv_avg';
%     h1 = plot(Pgeo(:,1),Pgeo(:,2),'k','linewidth',2); hold on; axis equal;
        
    % build gif
    figure(fig); frame = getframe(fig); 
    [A,map] = rgb2ind(frame2im(frame),256);
    if i == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.1);
    end
    delete([h1 h2]); 

    %% mesh AMG
    if msh_flg == 1
    % shift to center in mesh a bit
    AMG0(:,:,i) = [AMG0(:,1,i) + min(AMG0(:,1,i)), AMG0(:,2,i)];
    % write coordinates to *.dat file
    dlmwrite('airfoil.dat',[AMG0(:,:,i), zeros(N,1)],'delimiter','\t');

    %% Mesh airfoil using gmsh 2.10.1
    try 
    fprintf('calling dat2gmsh\n');
    !python dat2gmsh.py airfoil.dat
    %     !python dat2gmsh_v2.py airfoil.dat
    fprintf('calling gmsh\n');
    % may require sudo for encrypted hard drives
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

save('AMG_postproc.mat');