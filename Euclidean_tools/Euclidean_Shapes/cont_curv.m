function kappa = cont_curv(tN,th_spl,a_spl,th_scl,a_scl,Minv)
%spline evaluations
a = ppval(a_spl,tN);
th = ppval(th_spl,tN);

% spline derivatives
da  = a_scl*ppval(fnder(a_spl,1),tN); dda  = a_scl^2*ppval(fnder(a_spl,2),tN);
dth  = th_scl*ppval(fnder(th_spl,1),tN); ddth = th_scl^2*ppval(fnder(th_spl,2),tN);

% curvature
n1 = cos(th); n2 = sin(th);
dn1 = -sin(th); ddn1 = -cos(th);
dn2 = cos(th); ddn2 = -sin(th);
dc  = [da.*n1, da.*n2]*Minv' + ...
      [a.*dn1.*dth, a.*dn2.*dth]*Minv';
ddc = [dda.*n1, dda.*n2]*Minv' + ...
      2*[da.*dn1.*dth, da.*dn2.*dth]*Minv' + ...
      [a.*ddn1.*dth.^2, a.*ddn2.*dth.^2]*Minv' + ...
      [a.*dn1.*ddth, a.*dn2.*ddth]*Minv';
num = abs(dc(:,1).*ddc(:,2) - ddc(:,1).*dc(:,2));
den = ( dc(:,1).^2 + dc(:,2).^2 ).^(3/2);
kappa = num./den;