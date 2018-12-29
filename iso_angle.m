function Ptt0 = iso_angle(p,p0,Gt,U,Log)
% p: m x 1 point on the manifold
% p0: m x 1 point defining the central tangent space
% Gt: m x 1 tangent gradient vector at p
% U: m x (m-1) basis of the tangent space at p0
% Log: the inverse exponential map

% determine dimension
m = size(p,1);

% compute inverse exponentials of central and sample point
Vlog0 = -Log(p0,p); Vlog = Log(p,p0);

% compute the cosine and sine of the angle
cost = Vlog*Gt'/(sqrt(sum(Vlog.^2,2))*sqrt(sum(Gt.^2,2)));
sint = sqrt(1-cost^2);

% rotation in the subspace spanned by U
R = eye(m) - U(:,1)*U(:,1)' - U(:,2)*U(:,2)' + U*[cost -sint; sint cost]*U';
Ptt0 = sqrt(sum(Gt.^2,2))/sqrt(sum(Vlog0.^2,2))*Vlog0*R';