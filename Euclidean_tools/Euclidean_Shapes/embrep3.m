%  N: integer number of landmarks
% P0: matrix representation of shape landmarks

function [P_emb,nml,kappa,a_spl,th_spl,th0,a0,t,P_emb0,ind,tN] = embrep3(P0,N,smpl)
% number of points
n = size(P0,1);
%% affine transformations
% shift and scale to [-1,1] box
% pu = max(P0,[],1)'; pl = min(P0,[],1)';
% M = diag(2./(pu-pl)); b = -(M*pl + ones(2,1)); 

% shift and scale using landmark-affine standardization (Bryner, 2D affine and projective spaces)
% C_x = mean(P0,1)'; Sig_x = (P0 - repmat(C_x',n,1))'*(P0 - repmat(C_x',n,1));
% M = chol(Sig_x) \ eye(2); b = -M*C_x;
% % affine inverse (i.e., x0 = Minv*x + binv in [xl,xu]_m)
% Minv = chol(Sig_x); binv = C_x;

% shift and scale to center of mass
% pu = max(P0,[],1)'; pl = min(P0,[],1)';
% M = diag(2./(pu-pl));
% b = -M*mean(P0,1)';

% no transformation
M = eye(2); b = [0;0]; Minv = eye(2); ainv = b;

% transform points
P = (P0*M + repmat(b',n,1));

% repeat first point to close shape
P = [P;P(1,:)];

%% embedding
% compute discrete lengths
t0 = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 + ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
% compute unique angles using atan2
s0 = atan2(P(:,2),P(:,1)); 
s0 = unwrap(s0);

% circular embedding
nu_emb = [cos(s0) , sin(s0)];
% compute inner product
a0 = P(:,1).*nu_emb(:,1) + P(:,2).*nu_emb(:,2);
% take only unique points
[t0,ia] = unique(t0,'stable');
a0 = a0(ia); s0 = s0(ia);

% sort based on increasing discrete length
imgsort = sortrows([s0,a0,t0],3);
th0 = imgsort(:,1); a0 = imgsort(:,2); t0 = imgsort(:,3);

% scale length domain
lb = min(t0); ub = max(t0);
t = 1/(ub-lb)*t0;
% chain rule scaling
s_scl = 1/( (ub-lb) ); a_scl = 1/(ub-lb);

% build splines
% periodic spline of inner product
a_spl = csape(t,a0,'periodic');
% match first & second derivatives at endpoints of angular function
th_spl = csape(t,th0,'complete');

% find non-radially convex points
ds = ppval(th_spl,t);
ds = ds(2:end) - ds(1:end-1);
ind = (ds <= 0);

% evaluate shape at nominal landmarks
P_emb0 =  [ppval(a_spl,t).*cos(ppval(th_spl,t)),...
       ppval(a_spl,t).*sin(ppval(th_spl,t))];

% reevaluate shape at N landmarks using embedding representation
if strcmp(smpl,'curv')
    % compute continuous curvature approximation for refinement
    tcurv =linspace(0,1,10000)';
    da = a_scl*ppval(fnder(a_spl,1),tcurv); dda  = a_scl^2*ppval(fnder(a_spl,2),tcurv);
    ds = s_scl*ppval(fnder(th_spl,1),tcurv); dds = s_scl^2*ppval(fnder(th_spl,2),tcurv);
    num =  abs(-ppval(a_spl,tcurv).*dda.*ds + ppval(a_spl,tcurv).*da.*dds + 2*da.^2.*ds + ppval(a_spl,tcurv).^2.*ds.^3);
    den = ( da.^2 + (ppval(a_spl,tcurv).^2).*(ds.^2)).^(3/2);
    kappa = num./den;

    % refine domain using curvature based importance sampling
    ksum = cumsum(kappa); lb_k = min(ksum); ub_k = max(ksum);
    tN = pchip(1/(ub_k - lb_k)*ksum,tcurv,linspace(0,1,N)');

elseif strcmp(smpl,'uni')
    % uniform sampling for refinement
    tN = linspace(0,1,N)';
end

% compute continuous curvature at updated points 
% derivatives of splines
da = a_scl*ppval(fnder(a_spl,1),tN); dda  = a_scl^2*ppval(fnder(a_spl,2),tN);
ds = s_scl*ppval(fnder(th_spl,1),tN); dds = s_scl^2*ppval(fnder(th_spl,2),tN);
num =  abs(-ppval(a_spl,tN).*dda.*ds + ppval(a_spl,tN).*da.*dds + 2*da.^2.*ds + ppval(a_spl,tN).^2.*ds.^3);
den = ( da.^2 + ppval(a_spl,tN).^2.*ds.^2).^(3/2);
kappa = num./den;

% compute continuous normal vector approximation
nml  = a_scl*ppval(fnder(a_spl,1),tN).*([sin(ppval(th_spl,tN)), -cos(ppval(th_spl,tN))]) +...
      s_scl*ppval(fnder(th_spl,1),tN).*[ppval(a_spl,tN).*cos(ppval(th_spl,tN)), ppval(a_spl,tN).*sin(ppval(th_spl,tN))];
% compute unit normals
% nml = nml./sqrt(nml(:,1).^2 + nml(:,2).^2);

% re-evaluate shape at new N-landmarks
P_emb = [ppval(a_spl,tN).*cos(ppval(th_spl,tN)), ...
         ppval(a_spl,tN).*sin(ppval(th_spl,tN))];

%% transform back to original coordinates
% shape landmarks
P_emb = P_emb*Minv + repmat(ainv',N,1);
% plot continuous approximation of normals
nml = nml*[0 -1; 1 0]*Minv*[0 1; -1 0];