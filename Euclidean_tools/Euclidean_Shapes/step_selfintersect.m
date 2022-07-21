function [bool,p,q] = step_selfintersect(X)
% X is an n-by-2 discrete planar curve with n landmarks
% The last point must be duplicated as the first point for closed curves
% i.e., X(1,:) = X(n,:) is the "closure condition"

bool = false;

% infer number of landmarks
n = size(X,1);

% build diff operator
% Ds = -eye(n) + [zeros(1,n); [eye(n-1) zeros(n-1,1)]];
% diffX = X'*Ds;
% or use Matlab internal routine
diffX = diff(X)';
for p=1:n-2
    for q=(p+1):n-1
        A = [diffX(:,p) -diffX(:,q)];
        b = X(q,:)' - X(p,:)';
        invA = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))*[A(2,2) -A(1,2); -A(2,1) A(1,1)]; t = invA*b;
        if ( (10000*eps < t(1)) && (t(1) <1-10000*eps) ) && ( (10000*eps < t(2)) && (t(2) <1-10000*eps) )
            bool = true;
            return;
        end
    end
end