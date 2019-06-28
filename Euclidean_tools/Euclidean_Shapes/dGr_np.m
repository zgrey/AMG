function d = dGr_np(P1,P2)

% get dimensions (must be identical for P2)
p = size(P1,2);

% check if P1 and P2 are element of the Grassmannian
if abs(sum(diag(P1'*P1))/p - 1) >= 1e-8
    disp('WARNING: First input does not constitute an element of the Grassmannian')
elseif abs(sum(diag(P2'*P2))/p - 1) >= 1e-8
    disp('WARNING: Second input does not constitute an element of the Grassmannian')
end

% compute principal angles between subspaces
D = svd(P1'*P2,0);
theta = acos(D);

d = sqrt(sum(real(theta).^2));
