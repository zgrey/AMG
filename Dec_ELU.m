function Fout = Dec_ELU(X,W,b)

% ELU alpha term
alpha = 1;

% assign inputs
sig = X;

% depth
d = size(W,3);
for l=1:d
    if min(sum( abs(W(:,:,l)'), 1 )) ~= 0
%         sig = max( sig*W(:,:,l) + repmat(b(:,l)', size(sig,1), 1), zeros(size(sig,1),size(X,2)) + e);
        sig = max( sig*W(:,:,l)' + repmat(b(:,l)', size(sig,1), 1),...
                   alpha*(exp(sig*W(:,:,l)' + repmat(b(:,l)', size(sig,1), 1)) - 1));
    elseif min(sum( abs(W(:,:,l)), 1 )) == 0
        [~,ind] = max(sum( abs(W(:,:,l)'), 1), [], 2);
        [~,bi] = max(abs(b(:,l)));
%         sig = max( sig*W(:,ind,l) + repmat(b(bi,l), size(sig,1), 1), zeros(size(sig,1),1) + e);
        sig = max( sig*W(:,ind,l)' + repmat(b(bi,l), size(sig,1), 1),...
                   alpha*(exp(sig*W(:,ind,l)' + repmat(b(bi,l), size(sig,1), 1)) - 1));
    end
end

% return outputs
Fout = sig;