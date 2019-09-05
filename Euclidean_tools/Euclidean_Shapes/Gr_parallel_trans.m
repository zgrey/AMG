function H0 = Gr_parallel_trans(t,Y,H,Del)

% parallel translation of "Del" along geodesic in direction "H" a distance
% "t" starting from "Y"
% Edelman et al. (Theorem 2.4, pp. 321)
n = size(H,1);
[U,S,V] = svd(H,'econ');
H0 = ([Y*V U]*[-diag(sin(diag(S)*t) ); diag( cos(diag(S)*t) )]*U' + (eye(n) - U*U') )*Del;