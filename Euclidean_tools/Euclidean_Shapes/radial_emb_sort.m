% uni-variate radially convex embedding sort
function [ind,S,nu,S_rc,pp,s0,a0] = radial_emb_sort(s0,a0,N,lambda,CW,a,b)

% number of points
n = length(s0);

if CW == 1
    % clockwise sort
    imgsort = sortrows([-s0+1,a0],1);
else
    % counter-clockwise sort
    imgsort = sortrows([s0,a0],1);
end

s0 = imgsort(:,1); a0 = imgsort(:,2);

dS = repmat(s0,1,n) - repmat(s0',n,1);
dA = repmat(a0,1,n) - repmat(a0',n,1);

% penalized nearest-neighbor
R  = sqrt(dS.^2 + dA.^2) + lambda*dS;

R1S = repmat(s0,1,n);
[~,ISRT] = sort(R);
dSNRST = R1S(ISRT(2:end,:)) - repmat(s0',n-1,1);

ss = -ones(n,1); aa = -ones(n,1);
ss(1) = s0(1); aa(1) = a0(1); k = 1;
while k < n
    k = ISRT(find(dSNRST(:,k) > 0,1,'first')+1,k);
    
    ss(k) = s0(k);
    aa(k) = a0(k);
end
ind = find(ss == -1);
ss  = ss(ss ~= -1);
aa  = aa(aa ~= -1);
pp  = pchip([ss-1; ss; ss+1],[aa; aa; aa]);

if CW == 1
    % clockwise sort
    S =  [a*ppval(pp,s0).*cos((-s0+1)*2*pi), b*ppval(pp,s0).*sin((-s0+1)*2*pi)]; S = [S; S(1,:)];
else
    % CCW
    S =  [a*ppval(pp,s0).*cos(2*pi*(s0-1/2)), b*ppval(pp,s0).*sin(2*pi*(s0-1/2))]; S = [S; S(1,:)];
end

% compute discrete normal vectors
nu  = ppval(fnder(pp,1),linspace(0,1,N)').*([b*sin(2*pi*(linspace(0,1,N)'-1/2)), -a*cos(2*pi*(linspace(0,1,N)'-1/2))]) +...
     [b*2*pi*ppval(pp,linspace(0,1,N)').*cos(2*pi*(linspace(0,1,N)'-1/2)), a*2*pi*ppval(pp,linspace(0,1,N)').*sin(2*pi*(linspace(0,1,N)'-1/2))];
% compute discrete unit normals
nu = nu./sqrt(nu(:,1).^2 + nu(:,2).^2);

% re-evaluate the radially convex shape 
S_rc = [a*ppval(pp,linspace(0,1,N)').*cos((linspace(0,1,N)-1/2)*2*pi)', ...
        b*ppval(pp,linspace(0,1,N)').*sin((linspace(0,1,N)-1/2)*2*pi)'];