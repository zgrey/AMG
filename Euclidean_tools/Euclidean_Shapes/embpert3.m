function emb = embpert3(emb,del,A,B,type)

if strcmpi(type,'preserve')
    % volume preserving
    cvec = cos(emb.t*(1:length(A))*2*pi)*A;
    svec = sin(emb.t*(1:length(B))*2*pi)*B;
elseif strcmpi(type,'expand')
    cvec = cos(emb.t*(1:length(A))*2*pi)*A;
    svec = sin(emb.t*(1:length(B))*2*pi)*B;
    % volume expanding
    cvec = cvec + abs(min(cvec));
    svec = svec + abs(min(svec));
end

% perturb points
emb.pert.pts = emb.pts + del*(cvec + svec).*emb.nml;

tmp = embrep3(emb.pert.pts,emb.num,'uni');
emb.pert.nml  = tmp.nml;
emb.pert.curv = tmp.curv;
emb.pert.pts  = tmp.pts; 