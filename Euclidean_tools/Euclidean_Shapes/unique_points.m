
function [P,n,rep] = unique_points(P0)

%number of points
n0 = size(P0,1);
if n0 == 1
    P = P0; n = n0; rep =0;
else
    % sort for unique points
    D = squareform(pdist(P0)); rep = 0; ind = ones(1,n0);
    for i=1:n0
        ind = logical([ones(1,i),~(D(i,i+1:end) <= 1e-8)].*ind);
        rep = sum(D(i,i+1:end) <= 1e-8) + rep;   
    end
    P = P0(ind,:);
    % number of unique points
    n = size(P,1);
end