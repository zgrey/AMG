function V1 = schilds_ladder(p0,p1,V0,Exp,Log,N)
if norm(p0-p1) < 1e-10
    V1 = V0;
else
% discretize geodesic
t = linspace(0,1,N)'; v = Log(p0,p1); 
geo = Exp(t*norm(v),v,p0); 
% initialize
V1 = V0;

for i=1:N-1

    % sample geodesic discretization
    P0 = geo(i,:); P1 = geo(i+1,:);

    % Schild's ladder
    % scale rungs on ladder
    tau = t(i+1) - t(i);
    % geodesic from first point in direction
    P2  = Exp(1*tau,V1,P0);
    % "mid-point" between P2 and P1
    p   = Exp(0.5*tau,Log(P1,P2),P1);
    % vector mapping to twice the length through mid-point
    V1  = 1/tau*Log(P1,Exp(2*tau,Log(P0,p),P0));

end

end