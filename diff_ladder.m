function V1 = diff_ladder(p0,p1,V0,Exp,Log,N)
% discretize geodesic
t = linspace(0,1,N)'; v = Log(p0,p1);
geo = Exp(t*norm(v),v,p0); 
% initialize
V1 = V0;

for i=1:N-1

    % sample geodesic discretization
    P0 = geo(i,:); P1 = geo(i+1,:);

    % central-differencing ladder
    tau = 0.01;
    P2    = Exp(tau,V1,P0); P3    = Exp(-tau,V1,P0);
    Vlog2 = Log(P1,P2);     Vlog3 = Log(P1,P3);
    V1 = 1/(2*tau)*(Vlog2 - Vlog3);

end