function [H,B,U,S,V] = Gr_log(X,Y)

n = size(X,1);
p = size(X,2);

% Wang et al.
[U,S,V] = svd((Y - X*X'*Y) / (X'*Y),'econ');
H = U*diag(atan(diag(S)))*V';

% Rentmeesters et al., Gallivan et al., Shrivastava et al.
[W1,W2,~,C,~] = gsvd(X'*Y,Y - X*X'*Y,0);
A = W2(1:n-p,:)*diag(acos(diag(C)))*W1';
B = [zeros(p,p) -A';A zeros(n-p,n-p)];