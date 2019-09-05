function incr = mono_vec(F,varargin)

if ~isempty(varargin)
    tol = varargin{1};
else
    tol = 1e-6;
end

k = 1; incr = [1;zeros(length(F)-1,1)];
while k < length(F)
    if F(k) < F(k+1)
        k = k + 1;
        incr(k) = 1;
    else
        % skip until the next increasing value is found
        skp = find(F(k+1:end) - F(k) > tol,1);
        k = k + skp;
        incr(k) = 1;
    end
end
incr = logical(incr);
