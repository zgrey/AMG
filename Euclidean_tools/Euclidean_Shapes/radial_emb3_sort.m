% uni-variate radially convex embedding sort
function [ind,S,nu,S_emb,alpha,ss,s,a,t,kappa] = radial_emb3_sort(s0,a0,t0,N,smpl)

% sort based on increasing discrete length
imgsort = sortrows([s0,a0,t0],3);
s = imgsort(:,1); a = imgsort(:,2); t0 = imgsort(:,3);
lb = min(t0); ub = max(t0);

% scale length domain
t = 1/(ub-lb)*t0;

% build splines
% simple pchip's (not periodic)
% ss  = pchip(t,s);
% alpha  = pchip(t,a);
% periodic splines
ss = csape(t,s,'periodic');
alpha = csape(t,a,'periodic');
% match first & second derivatives at endpoints
% ss = csape(t,s,'second');
% alpha = csape(t,a,'second');
% match first derivative at endpoints
% ss = csape(t,s,'complete');
% alpha = csape(t,a,'complete');

% find non-radially convex points
ds = ppval(ss,t);
ds = ds(2:end) - ds(1:end-1);
ind = (ds <= 0);

% evaluate shape at nominal landmarks
S =  [ppval(alpha,t).*cos(2*pi*(ppval(ss,t)-1/2)),...
      ppval(alpha,t).*sin(2*pi*(ppval(ss,t)-1/2))];
% close shape
% S = [S; S(1,:)];

if strcmp(smpl,'curv')
% compute continuous curvature approximation for refinement
tcurv =linspace(0,1,10000)'; a_scl = 1/(ub-lb); s_scl = 1/(ub-lb);
da = a_scl*ppval(fnder(alpha,1),tcurv); dda  = a_scl^2*ppval(fnder(alpha,2),tcurv);
ds = s_scl*ppval(fnder(ss,1),tcurv); dds = s_scl^2*ppval(fnder(ss,2),tcurv);
num =  abs(-ppval(alpha,tcurv).*dda.*ds + ppval(alpha,tcurv).*da.*dds + 2*da.^2.*ds + ppval(alpha,tcurv).^2.*ds.^3);
den = ( da.^2 + (ppval(alpha,tcurv).^2).*(ds.^2)).^(3/2);
kappa = num./den;

% refine domain using curvature based importance sampling
ksum = cumsum(kappa); lb_k = min(ksum); ub_k = max(ksum);
tt = pchip(1/(ub_k - lb_k)*ksum,tcurv,linspace(0,1,N)');
elseif strcmp(smpl,'uni')
% uniform sampling
tt = linspace(0,1,N)';
elseif strcmp(smpl,'invcurv')
% refine domain using inverse curvature based importance sampling
ksum = cumsum(1./kappa); lb_k = min(ksum); ub_k = max(ksum);
tt = pchip(1/(ub_k - lb_k)*ksum,tt,linspace(0,1,N)');
end

% compute continuous curvature at updated points 
% chain rule scaling
s_scl = 1/( (ub-lb)*2*pi ); a_scl = 1/(ub-lb);
% derivatives of splines
da = a_scl*ppval(fnder(alpha,1),tt); dda  = a_scl^2*ppval(fnder(alpha,2),tt);
ds = s_scl*ppval(fnder(ss,1),tt); dds = s_scl^2*ppval(fnder(ss,2),tt);
num =  abs(-ppval(alpha,tt).*dda.*ds + ppval(alpha,tt).*da.*dds + 2*da.^2.*ds + ppval(alpha,tt).^2.*ds.^3);
% num = abs(ppval(alpha,tt).*ds.*(ppval(alpha,tt).^(ds.^2) - dda) + ppval(alpha,tt).*da.*dds + 2*(da.^2).*ds);
den = ( da.^2 + ppval(alpha,tt).^2.*ds.^2).^(3/2);
kappa = num./den;

% compute continuous normal vector approximation
nu  = ppval(fnder(alpha,1),tt).*([sin(2*pi*(ppval(ss,tt)-1/2)), -cos(2*pi*(ppval(ss,tt)-1/2))]) +...
     2*pi*ppval(fnder(ss,1),tt).*[ppval(alpha,tt).*cos(2*pi*(ppval(ss,tt)-1/2)), ppval(alpha,tt).*sin(2*pi*(ppval(ss,tt)-1/2))];
% compute unit normals
nu = nu./sqrt(nu(:,1).^2 + nu(:,2).^2);

% re-evaluate shape at new N-landmarks
S_emb = [ppval(alpha,tt).*cos((ppval(ss,tt)-1/2)*2*pi), ...
         ppval(alpha,tt).*sin((ppval(ss,tt)-1/2)*2*pi)];
% close shape
% S_emb = [S_emb; S_emb(1,:)]; 