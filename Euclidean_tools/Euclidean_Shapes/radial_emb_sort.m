% uni-variate radially convex embedding sort
function [ind,S,n,S_rc,pp] = radial_emb_sort(s0,a0,N,lambda,CW)

% number of points
l = length(s0);

if CW == 1
    % clockwise sort
    imgsort = sortrows([-s0+1,a0],1);
else
    % counter-clockwise sort
    imgsort = sortrows([s0,a0],1);
end

s0 = imgsort(:,1); a0 = imgsort(:,2);

dS = repmat(s0,1,l) - repmat(s0',l,1);
dA = repmat(a0,1,l) - repmat(a0',l,1);

% penalized nearest-neighbor
R  = sqrt(dS.^2 + dA.^2) + lambda*dS;

R1S = repmat(s0,1,l);
[~,ISRT] = sort(R);
dSNRST = R1S(ISRT(2:end,:)) - repmat(s0',l-1,1);

ss = -ones(l,1); aa = -ones(l,1);
ss(1) = s0(1); aa(1) = a0(1); k = 1;
while k < l
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
    S =  [ppval(pp,s0).*cos((-s0+1)*2*pi), ppval(pp,s0).*sin((-s0+1)*2*pi)]; S = [S; S(1,:)];
else
    % CCW
    S =  [ppval(pp,s0).*cos(s0*2*pi), ppval(pp,s0).*sin(s0*2*pi)]; S = [S; S(1,:)];
end

% compute discrete normal vectors
% 90 degree rotation
G = [0 -1;1 0];
% compute unit normals
n  = ppval(fnder(pp,1),s0).*([cos(2*pi*s0),sin(2*pi*s0)]*G') - [2*pi*ppval(pp,s0).*cos(2*pi*s0), 2*pi*ppval(pp,s0).*sin(2*pi*s0)];
n = -n./sqrt(n(:,1).^2 + n(:,2).^2);

% re-evaluate the radially convex shape 
S_rc = [ppval(pp,linspace(0,1,N)').*cos(linspace(0,1,N)*2*pi)', ...
        ppval(pp,linspace(0,1,N)').*sin(linspace(0,1,N)*2*pi)'];