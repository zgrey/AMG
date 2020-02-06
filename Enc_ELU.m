function Fout = Enc_ELU(X,W,b)

% "leaky" term
e = 0;

% ELU alpha term
alpha = 1;

% assign inputs
sig = X;

% depth
d = size(W,3);
for l=1:d
    % identify non-zero entries
    if min(sum( abs(W(:,:,l)), 1 )) ~= 0
        % non-smooth "rectified linear unit"
%         sig = max( sig*W(:,:,l) + repmat(b(:,l)', size(sig,1), 1), zeros(size(sig,1),size(X,2)) + e);
        % smooth exponential "linear unit"
        sig = max( sig*W(:,:,l) + repmat(b(:,l)', size(sig,1), 1),...
                   alpha*(exp(sig*W(:,:,l) + repmat(b(:,l)', size(sig,1), 1)) - 1));
    elseif min(sum( abs(W(:,:,l)), 1 )) == 0
        % obtain non-zero indices
        [~,ind] = max(sum( abs(W(:,:,l)), 1), [], 2);
        [~,bi] = max(abs(b(:,l)));
        % non-smooth "rectified linear unit"
%         sig = max( sig*W(:,ind,l) + repmat(b(bi,l), size(sig,1), 1), zeros(size(sig,1),1) + e);
        % smooth exponential "linear unit"
        sig = max( sig*W(:,ind,l) + repmat(b(bi,l), size(sig,1), 1),...
                   alpha*(exp(sig*W(:,ind,l) + repmat(b(bi,l), size(sig,1), 1)) - 1));
    end
end

% return outputs
Fout = sig;