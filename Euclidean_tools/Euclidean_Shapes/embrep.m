% Compute circular embedding pchip representation:
%
% Inputs:
% C - R^(n x 2) matrix representing a subset of particles from a shape (no
%     duplicate particles/landmarks)
% N - number of particles/points from the embedding representation
% t - scalar perturbation parameter
% A - tall vector of cosine coefficients
% B - tall vector of sine coefficients
% NOTE: the points C are assumed to be feasible for circular embedding 
% (i.e., near circle)
%
% Outputs:
%   vs - R^n vector containing the embedding function values at the n particles
%    s - R^n vector containing the domain values for the function at the n particles
% nemb - R^(n x 2) matrix containing scaled unit normal vectors at the n particles
%   CC - R^(N x 2) new coordinates from the approximated embedding function
%    k - R^N curvature values calculated at the approximated N coordinates
%   ss - R^N vector of new uniformly distributed domain values
%   vv - R^N vector of approximated embedding function values at entries of ss
%   nn - R^(N-1) vector of unit normals

% 
function [vs,s,nemb,CC,k,ss,vv,nn,sarc] = embrep(C,N,t,A,B)

%% determine nominal particle discrete embedding
% duplicate first landmark from data
C = [C; C(1,:)];
theta = atan2(C(1:end-1,2),C(1:end-1,1)); negi = find(theta < 0);
theta(negi) = theta(negi) + 2*pi;
% unit normals of circle
nc = [cos(theta), sin(theta)];
% recover each landmark
nemb = sqrt(C(1:end-1,1).^2+C(1:end-1,2).^2).*nc;
% discrete embedding representation
vs = sum(C(1:end-1,:).*nc,2); cssort = sortrows([theta,vs],1);
vs = cssort(:,2); s = cssort(:,1)/(2*pi);

%% domain distribution
% curvature based sampling
[CDF] = ecdf(tesscurv(C)); ss = pchip(CDF,[0;s],linspace(0,1,N)');
% if isempty(find(ss==0, 1)), ss = [0;ss]; N=N+1; end
% if isempty(find(ss==1, 1)), ss = [ss;1]; N=N+1; end

% pchip of cummulative sum for non-uniform distribution of domain samples
% vv = pchip([s-1;s;1+s],[vs;vs;vs],[0;s;1]); csum = cumsum(vv);
% ssinv = pchip((csum-min(csum))/(max(csum)-min(csum)),[0;s;1],linspace(0,1,N)'); ss = sort(ssinv);

% unifrom distribution of 
% ss = linspace(0,1,N)';

% empirical cdf of alpha
% [CDF,vcdf] = ecdf(vs); ssinv = pchip(CDF,vcdf,linspace(0,1,N)'); ss = ssinv/(max(ssinv)-min(ssinv)) -min(ssinv);

%% refinement of particles/landmarks using periodic approximation 
% 3-period pchip
vv = pchip([s-1;s;1+s],[vs;vs;vs],ss); vp = pchip([s-1;s;1+s],[vs;vs;vs]);

% periodic spline
% vp = cscvn([s';vs']); vvp = fnval(vp,pchip(s,vp.breaks,ss));
% vv = vvp(2,:)'; ss = vvp(1,:)';

% refined particles/landmarks
CC = [vv.*cos(ss*2*pi), vv.*sin(ss*2*pi)];
% tangent vectors
tau  = (CC(2:end,:) - CC(1:end-1,:));

% discrete normal vectors
% nn(1,:) = 1/2*(CC(2,:)-CC(N-1,:))*[0 -1;1 0];
% nn(2:N-1,:) = 1/2*(CC(3:N,:)-CC(1:N-2,:))*[0 -1;1 0]; nn = nn./sqrt(sum(nn.*nn,2));
% continuous normal vectors
nn = ppval(fnder(vp,1),ss).*([cos(2*pi*ss),sin(2*pi*ss)]*[0 -1;1 0]') - ...
    [2*pi*ppval(vp,ss).*cos(2*pi*ss), 2*pi*ppval(vp,ss).*sin(2*pi*ss)];
nn = -nn./sqrt(nn(:,1).^2 + nn(:,2).^2);

%% discrete normal perturbation to shape
sarc = [0;cumsum(sqrt(tau(:,1).^2 + tau(:,2).^2))]/sum(sqrt(tau(:,1).^2 + tau(:,2).^2));
% area preserving
cvec = cos(sarc*(1:length(A))*2*pi)*A;
svec = sin(sarc*(1:length(B))*2*pi)*B;
% area expanding
cvec = cvec + abs(min(cvec));
svec = svec + abs(min(svec));
% smooth normal field
CC = CC + t*(cvec + svec).*nn; CC(end,:) = CC(1,:);
% discrete normal vectors
nn(1,:) = 1/2*(CC(2,:)-CC(N-1,:))*[0 -1;1 0];
nn(2:N-1,:) = 1/2*(CC(3:N,:)-CC(1:N-2,:))*[0 -1;1 0]; nn = nn./sqrt(sum(nn.*nn,2));
% discrete curvature computation
% Ref: Differential Geometry of Curves and Surfaces, do Carmo (Ch. 1)
% Ref: Differential Geometry of Curves and Surfaces, Toponogov
% discrete tesselation curvature of N approximated points
ktess = tesscurv(CC);

% p-discrete curvature
% Ref: Parabola-Based Discrete Curvature Estimation, Kim, H., Rossignac, J.
CC0 = [CC, zeros(N,1)];
CCfwd = [CC0(2:end,:) ; CC0(2,:)]; CCbwd = [CC0(end-1,:) ; CC0(1:end-1,:)];
area = cross(4*(CCbwd-2*CC0+CCfwd),CCfwd-CCbwd); step = CCfwd - CCbwd;
kp = abs(area(:,3)) ./ ((step(:,1).^2 + step(:,2).^2).^(3/2));

% Pick curvature
k = ktess;
% k = kp;

%% emedding perturbation
% sarc = [0;cumsum(sqrt(tau(:,1).^2 + tau(:,2).^2))]/sum(sqrt(tau(:,1).^2 + tau(:,2).^2));
% % volume preserving
% cvec = cos(sarc*(1:length(A))*2*pi)*A;
% svec = sin(sarc*(1:length(B))*2*pi)*B;
% % volume expanding
% cvec = cvec + (sum(A))*ones(length(s),1);
% svec = svec + (sum(B))*ones(length(s),1);
% vv = vv + t*(cvec + svec);
% vp = pchip([ss-1;ss;1+ss],[vv;vv;vv]); vv = ppval(vp,ss);
% % landmarks
% CC = [vv.*cos(ss*2*pi), vv.*sin(ss*2*pi)];
% % continuous normal vectors
% nn = ppval(fnder(vp,1),ss).*([cos(2*pi*ss),sin(2*pi*ss)]*[0 -1;1 0]') - ...
%     [2*pi*ppval(vp,ss).*cos(2*pi*ss), 2*pi*ppval(vp,ss).*sin(2*pi*ss)];
% nn = -nn./sqrt(nn(:,1).^2 + nn(:,2).^2);
% embedding representation curvature
% n1   = cos(2*pi*ss); n2 = sin(2*pi*ss);
% dn1  = -sin(2*pi*ss)*2*pi; dn2 = cos(2*pi*ss)*2*pi;
% ddn1 = -cos(2*pi*ss)*4*pi^2; ddn2 = -sin(2*pi*ss)*4*pi^2;
% dv   = ppval(fnder(vp,1),ss); ddv  = ppval(fnder(vp,2),ss);
% kemb = abs( (dv.*n1 + vv.*dn1).*(ddv.*n2 + 2*dv.*dn2 + vv.*ddn2) - ...
%             (dv.*n2 + vv.*dn2).*(ddv.*n1 + 2*dv.*dn1 + vv.*ddn1) ) ./ ...
%             ( (dv.*n1 + vv.*dn1).^2 + (dv.*n2 + vv.*dn2).^2 ).^(3/2);
% k = kemb;
