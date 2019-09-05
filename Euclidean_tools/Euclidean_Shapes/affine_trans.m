function [P_TF,M,b,Minv] = affine_trans(P,type)

% number of landmarks
n = size(P,1);

if strcmpi(type,'cube')
    % shift and scale to [-1,1] box
    pu = max(P,[],1)'; pl = min(P,[],1)';
    % affine transformation
    M = diag(2./(pu-pl)); Minv = diag((pu-pl)/2);
    b = -(M*pl + ones(2,1)); 
    
elseif strcmpi(type,'LA')
    % shift and scale using landmark-affine standardization (Bryner, 2D affine and projective spaces)
    C_x = mean(P,1)';
    [U,D] = svd(1/sqrt(n-1)*(P - repmat(C_x',n,1))',0); D = D(1:2,1:2);
    M = diag(1./diag(D))*U'; Minv = U*D;

%     Sig_x = 1/(n-1)*(P - repmat(C_x',n,1))'*(P - repmat(C_x',n,1)); [U,D] = svd(Sig_x);
%     Sig_x = 1/(n-1)*(P - repmat(C_x',n,1))'*(P - repmat(C_x',n,1)); [U,D] = eigs(Sig_x);
%     [U,D] = svd(1/sqrt(n-1)*(P - repmat(C_x',n,1))',0); D = D(1:2,1:2).^2;
    % affine transformation
%     M = diag(1./sqrt(diag(D)))*U'; Minv = U*sqrt(D);
    b = -M*C_x;

elseif strcmpi(type,'CM')
    % shift and scale to center of mass
    pu = max(P,[],1)'; pl = min(P,[],1)';
    % affine transformation
    M = diag(2./(pu-pl)); Minv = diag((pu-pl)/2);
    b = -M*mean(P,1)';
    
elseif strcmpi(type,'none')
    % no transformation
    M = eye(2); Minv = eye(2);
    b = [0;0]; 
    
end

P_TF = P*M' + repmat(b',n,1);
