% uni-variate radially convex embedding sort
function [ind,S,nu,S_emb,pp,ss,s,a,t,kappa] = radial_emb3_sort(s0,a0,t0,N)

% sort based on increasing discrete length
imgsort = sortrows([s0,a0,t0],3);
s = imgsort(:,1); a = imgsort(:,2); t0 = imgsort(:,3);

% scale length domain (without extrapolation)
lb = min(t0); ub = max(t0);
t = 1/(ub-lb)*t0;

% build splines
% simple pchip's (not periodic)
% ss  = pchip(t,s);
% pp  = pchip(t,a);
% periodic splines
% ss = csape(t,s,'periodic');
% pp = csape(t,a,'periodic');
% match first & second derivatives at endpoints
ss = csape(t,s,'second');
pp = csape(t,a,'second');
% match first derivative at endpoints
% ss = csape(t,s,'complete');
% pp = csape(t,a,'complete');

% find non-radially convex points
ds = ppval(ss,t);
ds = ds(2:end) - ds(1:end-1);
ind = (ds <= 0);

% evaluate shape at nominal landmarks
S =  [ppval(pp,t).*cos(2*pi*(ppval(ss,t)-1/2)),...
      ppval(pp,t).*sin(2*pi*(ppval(ss,t)-1/2))];
% close shape
% S = [S; S(1,:)];

% compute continuous curvature approximation for refinement
tt =linspace(0,1,10000)';
n1   = cos(2*pi*(ppval(ss,tt)-1/2)); n2 = sin(2*pi*(ppval(ss,tt)-1/2));
dn1  = -2*pi*sin(2*pi*(ppval(ss,tt)-1/2)).*ppval(fnder(ss,1),tt); dn2 = 2*pi*cos(2*pi*(ppval(ss,tt)-1/2)).*ppval(fnder(ss,1),tt);
ddn1 = -4*pi^2*cos(2*pi*(ppval(ss,tt)-1/2)).*ppval(fnder(ss,2),tt); ddn2 = -4*pi^2*sin(2*pi*(ppval(ss,tt)-1/2)).*ppval(fnder(ss,2),tt);
dv   = ppval(fnder(pp,1),tt); ddv  = ppval(fnder(pp,2),tt);
kappa = abs( (dv.*n1 + ppval(pp,tt).*dn1).*(ddv.*n2 + 2*dv.*dn2 + ppval(pp,tt).*ddn2) - ...
            (dv.*n2 + ppval(pp,tt).*dn2).*(ddv.*n1 + 2*dv.*dn1 + ppval(pp,tt).*ddn1) ) ./ ...
            ( (dv.*n1 + ppval(pp,tt).*dn1).^2 + (dv.*n2 + ppval(pp,tt).*dn2).^2 ).^(3/2);

% refine domain using curvature based importance sampling
ksum = cumsum(kappa); lb_k = min(ksum); ub_k = max(ksum);
tt = pchip(1/(ub_k - lb_k)*ksum,tt,linspace(0,1,N)');
% uniform sampling
% tt = linspace(0,1,N)';

% compute continuous curvature
n1   = cos(2*pi*(ppval(ss,tt)-1/2)); n2 = sin(2*pi*(ppval(ss,tt)-1/2));
dn1  = -2*pi*sin(2*pi*(ppval(ss,tt)-1/2)).*ppval(fnder(ss,1),tt); dn2 = 2*pi*cos(2*pi*(ppval(ss,tt)-1/2)).*ppval(fnder(ss,1),tt);
ddn1 = -4*pi^2*cos(2*pi*(ppval(ss,tt)-1/2)).*ppval(fnder(ss,2),tt); ddn2 = -4*pi^2*sin(2*pi*(ppval(ss,tt)-1/2)).*ppval(fnder(ss,2),tt);
dv   = ppval(fnder(pp,1),tt); ddv  = ppval(fnder(pp,2),tt);
kappa = abs( (dv.*n1 + ppval(pp,tt).*dn1).*(ddv.*n2 + 2*dv.*dn2 + ppval(pp,tt).*ddn2) - ...
            (dv.*n2 + ppval(pp,tt).*dn2).*(ddv.*n1 + 2*dv.*dn1 + ppval(pp,tt).*ddn1) ) ./ ...
            ( (dv.*n1 + ppval(pp,tt).*dn1).^2 + (dv.*n2 + ppval(pp,tt).*dn2).^2 ).^(3/2);

% compute continuous normal vector approximation
nu  = ppval(fnder(pp,1),tt).*([sin(2*pi*(ppval(ss,tt)-1/2)), -cos(2*pi*(ppval(ss,tt)-1/2))]) +...
     2*pi*ppval(fnder(ss,1),tt).*[ppval(pp,tt).*cos(2*pi*(ppval(ss,tt)-1/2)), ppval(pp,tt).*sin(2*pi*(ppval(ss,tt)-1/2))];
% compute unit normals
nu = nu./sqrt(nu(:,1).^2 + nu(:,2).^2);

% re-evaluate shape at new N-landmarks
S_emb = [ppval(pp,tt).*cos((ppval(ss,tt)-1/2)*2*pi), ...
         ppval(pp,tt).*sin((ppval(ss,tt)-1/2)*2*pi)];
% close shape
% S_emb = [S_emb; S_emb(1,:)]; 