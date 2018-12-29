function V1 = diff_ladder(p0,p1,V0,Exp,Log,N)
if norm(p0-p1) < 1e-10
    V1 = V0;
else
% discretize geodesic
t = linspace(0,1,N)'; v = Log(p0,p1); 
if N ~= 1, dt = 1/(N-1); end
geo = Exp(t*norm(v),v,p0); 
% initialize
V1 = V0;

for i=1:N-1

    % sample geodesic discretization
    P0 = geo(i,:); P1 = geo(i+1,:);

    % central-differencing ladder
    tau = 1e-2;
    P2    = Exp(tau,V1,P0); P3    = Exp(-tau,V1,P0);
    Vlog2 = Log(P1,P2);     Vlog3 = Log(P1,P3);
    V1 = 1/(2*tau)*(Vlog2 - Vlog3);

end

end