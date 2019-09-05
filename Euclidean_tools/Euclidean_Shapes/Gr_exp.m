function Y = Gr_exp(t,X,H)

[U,S,V] = svd(H,'econ');

% Consistent between Edelman et al. and Wang et al. (seems to contrast
% Rentmeesters et al., looks like a typo bullet 2., pp. 3840)
% Y = X*V*diag(cos(t*diag(S)))*V' + U*diag(sin(t*diag(S)))*V';
Y = [X*V U]*[ diag(cos(t*diag(S))); diag(sin(t*diag(S))) ]*V';