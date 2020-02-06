function Fout = NN_ELU(X,Wvec,bvec,strc,alpha)

% initialize inputs
sig = X;

% depth
d = length(strc)-1; shft_W = 1; shft_b = 1;
for l=1:d
    % build affine transformation
    W = reshape(Wvec(shft_W:shft_W -1 + strc(l)*strc(l+1)),strc(l+1),strc(l));
    b = bvec(shft_b:shft_b - 1 + strc(l+1));
    % smooth exponential "linear unit"
    sig = max( sig*W' + repmat(b', size(sig,1), 1),...
               alpha*(exp(sig*W' + repmat(b', size(sig,1), 1)) - 1));
    % parameter indexing
    shft_W  = strc(l)*strc(l+1) + 1;
    shft_b = strc(l+1) + 1;
end

% return outputs
Fout = sig;