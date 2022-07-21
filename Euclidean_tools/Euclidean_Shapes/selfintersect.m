function [bool,Xint,ind,tp,tq] = selfintersect(X)
% X is an n-by-2 discrete planar curve with n landmarks
% The last point must be duplicated as the first point for closed curves
% i.e., X(1,:) = X(n,:) is the "closure condition"

bool = false;
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

% infer number of landmarks
n = size(X,1);

% build diff operator
% Ds = -eye(n) + [zeros(1,n); [eye(n-1) zeros(n-1,1)]];
% diffX = X'*Ds;
% or use Matlab internal routine
diffX = diff(X)';

ind = zeros(n,n);
tp = zeros(n,1); tq = zeros(n,1); Xint = zeros(n,2);
for p=1:n-2
    for q=(p+1):n-1
        A = [diffX(:,p) -diffX(:,q)];
        b = X(q,:)' - X(p,:)';
        invA = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))*[A(2,2) -A(1,2); -A(2,1) A(1,1)]; t = invA*b;
        if ( (10000*eps < t(1)) && (t(1) <1-10000*eps) ) && ( (10000*eps < t(2)) && (t(2) <1-10000*eps) )
            ind(p,q) = 1;
            tp(p) = t(1);
            tq(q) = t(2);
            Xint(p,:) = t(1)*X(p+1,:) + (1-t(1))*X(p,:);
        end
    end
end
Xint = Xint(tp ~= 0,:);

if sum(sum(ind)) > 0
    bool = true;
end

warning('on','MATLAB:singularMatrix');
warning('on','MATLAB:nearlySingularMatrix');