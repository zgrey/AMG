% Model Adjoint Sensitivities of QiQi airfoil database
% Zach Grey - 08/26/2019
clc; close all; clearvars;

addpath ~/AMG/Euclidean_tools/Euclidean_Shapes/

% number of resampled points
N = 10000; plt_flg = 0;

%% read nominal points from file
dirstruct  = dir('~/airfoil_coordinate_database/coordinates/');
Naf = size(dirstruct,1)-2; i_fail = zeros(Naf,1); mesh_fail = i_fail;

if plt_flg == 1, fig = figure('units','normalized','outerposition',[0 0 1 1]); end
for i=1:Naf
    % obtain nominal landmarks
    try 
        P0 = importairfoil([dirstruct(2+i).folder,'/',dirstruct(2+i).name]);
    catch
        disp(['Could not read: ',dirstruct(2+i).name])
        i_fail(i) = 1;
    end
    % number of points
    n = size(P0,1);


    %% compute embedding representation
    emb = embrep3(P0,N,'uni');

    %% Visualize air   foils
    if plt_flg == 1
        % plot continuous approximation of shape
        subplot(1,2,1), h1 = plot(emb.pts(:,1),emb.pts(:,2),'linewidth',2); axis equal; hold on;
        subplot(1,2,1), h2 = quiver(emb.pts(:,1),emb.pts(:,2),emb.nml(:,1),emb.nml(:,2),'color',h1.Color);
        % plot original airfoil
        subplot(1,2,1), h3 = scatter(P0(:,1),P0(:,2),15,'k');
        fig.CurrentAxes.Visible = 'off';
        % plot transformed coordinates
        subplot(1,2,2), h4 = plot(emb.TF.pts(:,1),emb.TF.pts(:,2),'linewidth',2); hold on;
        TF = emb.pts*emb.TF.M' + repmat(emb.TF.b',N,1);
        subplot(1,2,2), h7 = plot(TF(:,1),TF(:,2),'--','linewidth',2);
        subplot(1,2,2), h5 = scatter(emb.TF.pts(:,1),emb.TF.pts(:,2),15,'k');

        h6 = annotation('textbox', [0 0.9 1 0.1], ...
            'String', ['Directory index i=',num2str(i),', ',dirstruct(2+i).name], ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');

        % print figure
        print(['./figs_airfoils/',num2str(i),'_',dirstruct(2+i).name,'.png'],'-dpng')
        % clear plot
        delete([h1;h2;h3;h4;h5;h6;h7]); reset(gca);
        clc; disp([num2str(i),'/',num2str(Naf),' Shapes Complete']);
    end
    
    
    % write coordinates to *.dat file
    dlmwrite('airfoil.dat',[emb.pts, zeros(size(emb.pts,1),1)],'delimiter','\t');

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
    % save mesh for visualization
    copyfile('./airfoil.mesh',['./gmsh_meshes/airfoil_',num2str(i),'_',dirstruct(2+i).name,'.mesh']);

    %% Convert Mesh
    tic;
    fprintf('mesh %i of %i pre-processing/conversion to .su2 format...\n',i,Naf);
    % Run DARPA EQUiPS SEQUOIA team codes:
    meshGMF = ReadGMF('airfoil.mesh');
    meshGMF = meshPrepro(meshGMF);
    meshSU2 = convertGMFtoSU2(meshGMF);
    WriteSU2_airfoil(meshSU2, ['./meshes/airfoil_',num2str(i),'.su2']);

    tocp = toc;
    % print time stats
    disp(['Finished in... ',num2str(toc),' sec | Remaining time ~',num2str((Naf-i)*(toc+tocp)/120),' min']);
    catch
        mesh_fail(i) = 1;
    end
end
