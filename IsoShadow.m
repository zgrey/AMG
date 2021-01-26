
function ISO = IsoShadow(Y,F,Nlvl)
% A function for visualizing 3D level sets of a 3D ridge function
% Y:    N by 3 matrix where N are the number of samples and 3 is the number of
%       active parameters
% F:    N by 1 matrix function evaluations paired with rows of Y
% Nlvl: The number of level sets over the observed range of F

% start by partitioning the range of F
maxF = max(F);
minF = min(F);
% partition level sets
if Nlvl == 1
    lvlsets = mean(F);
else
    lvlsets = linspace(0.9*minF,0.9*maxF,Nlvl);
end

% cube the data for isosurface
Nshp = floor(size(Y,1)^(1/3)); YY = zeros(Nshp,Nshp,Nshp,3);

for i=1:3, YY(:,:,:,i) = reshape(Y(1:Nshp^3,i),Nshp,Nshp,Nshp); end

FV = cell(Nlvl,1);
for i=1:Nlvl
    fprintf('Building isosurface %i of %i...',i,Nlvl);
    FV{i} = isosurface(YY(:,:,:,1), YY(:,:,:,2), YY(:,:,:,3),reshape(F(1:Nshp^3),Nshp,Nshp,Nshp),lvlsets(i));
    Ptch = patch(FV{i},'FaceVertexCData',lvlsets(i)*ones(size(FV{i}.vertices,1),1),'FaceColor','flat'); 
    colorbar; grid on;
    Ptch.FaceAlpha = 0.5; Ptch.EdgeColor = 'none';
    clc;
end

ISO = FV;