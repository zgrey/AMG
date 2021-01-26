% Iso-Shadow plot demo (an alternative is slice() )
clc; close all;

%% Things to modify
% The number of original parameters (ambient dimension)
m = 10;
% ridge function type (lin or quad)
% type = 'lin';
type = 'quad';
% number of level-sets to plot
Nlvl = 20;
    
%% ridge functions
if strcmpi(type,'lin')
    % 3d linear ridge function
    A = 2*rand(m,3) -1;
    % active directions
    W = [A(:,1)/norm(A(:,1)), A(:,2)/norm(A(:,2)), A(:,3)/norm(A(:,3))];
    % function handle
    f = @(X) sum(X*W,2);
elseif strcmpi(type,'quad')
    % 3d nonlinear ridge function 
    A = 2*rand(m,3) -1;
    % active directions
    W = [A(:,1)/norm(A(:,1)), A(:,2)/norm(A(:,2)), A(:,3)/norm(A(:,3))];
    % function handle
    f = @(X) sum((X*[W zeros(m,m-3)]).*X,2);
end

%% Random domain
% some uniform [-1,1]^m random samples
N = 10000;
X = 2*rand(N,m) - 1;
% evaluate and project random data to 3d active subspace
F = f(X); Y = X*W;
% plot level sets
figure; IsoShadow(Y,F,Nlvl); view(-37.5,30);
title 'Random data'

%% Structured domain
% alternatively, a structured grid on R^3
Nshp = 100;
[Y1,Y2,Y3] = meshgrid(linspace(-1,1,Nshp),...
                      linspace(-1,1,Nshp),...
                      linspace(-1,1,Nshp));
Y = [reshape(Y1,Nshp^3,1), reshape(Y2,Nshp^3,1), reshape(Y3,Nshp^3,1)];
% evaluate along ridge (or 3d surrogate)
X = Y*W'; F = f(X);

% start by partitioning the range of F
maxF = max(F);
minF = min(F);
% assign appropriate grid manually
YY(:,:,:,1) = Y1;
YY(:,:,:,2) = Y2;
YY(:,:,:,3) = Y3;
% build uniform level-set values at 90% of min and max
lvlsets = linspace(0.9*minF,0.9*maxF,Nlvl);

FV = cell(Nlvl,1); figure; 
for i=1:Nlvl
    fprintf('Building isosurface %i of %i...',i,Nlvl);
    FV{i} = isosurface(YY(:,:,:,1), YY(:,:,:,2), YY(:,:,:,3),reshape(F(1:Nshp^3),Nshp,Nshp,Nshp),lvlsets(i));
    Ptch = patch(FV{i},'FaceVertexCData',lvlsets(i)*ones(size(FV{i}.vertices,1),1),'FaceColor','flat'); 
    colorbar; grid on;
    Ptch.FaceAlpha = 0.5; Ptch.EdgeColor = 'none';
    clc;
end
title 'Structured data'; view(-37.5,30);