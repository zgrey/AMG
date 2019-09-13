function emb = embrep3(P0,N,smpl,varargin)
%   P0: matrix representation of shape landmarks
%    N: integer number of landmarks
% smpl: sampling type ('uni' or 'curv')

% set default options
opts = struct('ThetaSpline','complete',...
              'AffineTrans','LA',...
              'M',[],...
              'b',[],...
              'Minv',[]);
          
if ~isempty(varargin)
    opts_custom = varargin{1};
    all_fields = fieldnames(opts);
    for field = fieldnames(opts_custom)'
        tf_fieldmatch = strcmpi(all_fields, field);
        if any(tf_fieldmatch)
            opts = setfield(opts, all_fields{tf_fieldmatch}, getfield(opts_custom, field{:}));
        else
            error(['ERROR: Unexpected field name -> ' field])
        end
    end
end

%% Affine transformations
% nominal amount of provided landmark data
emb.nom.num = size(P0,1);

% find unique points (very slow for large number of landmarks, currently limited by squareform)
[P,n] = unique_points(P0); emb.nom.unq_num = n; emb.nom.num = size(P0,1);
% SPEED-UPS FOR CONVERGENCE STUDY:
% if data has a known repeated last point (remove it for the transformation)
% P = P0(1:end-1,:); n = size(P0,1)-1;  emb.nom.unq_num = size(P0,1)-1;
% if data does not have a repeated last point (repeat the first for periodicity, line 67)
% P = P0; n = size(P0,1); emb.nom.unq_num = size(P0,1);

% compute affine transformations
if strcmpi(opts.AffineTrans,'custom')
    M = opts.M; b = opts.b; Minv = opts.Minv;
else
    [~,M,b,Minv] = affine_trans(P,opts.AffineTrans);
end

% transform points
P = P*M' + repmat(b',n,1);
P0 = P0*M' + repmat(b',emb.nom.num,1);
% repeat first unique point to close shape (periodicity requirement for spline)
P = [P;P(1,:)];

% save transformation
emb.TF.M = M; emb.TF.b = b; emb.TF.Minv = Minv; emb.TF.pts = P;

%% Embedding representation
% compute discrete lengths of non-repeated points
t  = cumsum([0; sqrt(( P(2:end,1) - P(1:end-1,1) ).^2 + ( P(2:end,2) - P(1:end-1,2) ).^2)],1);
t0 = cumsum([0; sqrt(( P0(2:end,1) - P0(1:end-1,1) ).^2 + ( P0(2:end,2) - P0(1:end-1,2) ).^2)],1);
% compute unique angles using atan2 of original landmarks
th = unwrap(atan2(P(:,2),P(:,1))); 
% circular embedding
circ_nml = [cos(th) , sin(th)];
% compute inner product
alph = P(:,1).*circ_nml(:,1) + P(:,2).*circ_nml(:,2);

% save transformed landmark angle and inner product
emb.TF.th = th; emb.TF.alph = alph;

% scale length domain
ub = max(t);
emb.TF.t = t/ub; emb.nom.t = t0/ub;
% chain rule scaling
emb.th.scl = 1/ub; emb.alph.scl = 1/ub;

% build splines
if strcmp(opts.ThetaSpline,'complete')
    % match endslopes at endpoints of angular function
    emb.th.spl = csape(emb.TF.t,emb.TF.th,'complete');
elseif strcmp(opts.ThetaSpline,'pchip')
    % respects strictly monotonic data
    emb.th.spl = pchip(emb.TF.t,emb.TF.th);
end
% periodic spline of inner product
emb.alph.spl = csape(emb.TF.t,emb.TF.alph,'periodic');

%% Sampling
% reevaluate shape at N landmarks using embedding representation
if strcmp(smpl,'curv')
    % compute continuous curvature approximation for refinement
    tcurv = linspace(0,1,10000)';
    emb.curv = cont_curv(tcurv,emb.th.spl,emb.alph.spl,emb.th.scl,emb.alph.scl,Minv);

    % refine domain using curvature-based importance sampling
    ksum = cumsum(emb.curv); lb_k = min(ksum); ub_k = max(ksum);
    tN = pchip(1/(ub_k - lb_k)*ksum,tcurv,linspace(0,1,N)');

elseif strcmp(smpl,'uni')
    tN = linspace(0,1,10000)';
    % evaluate shape at new N-landmarks
    tmp_pts = [ppval(emb.alph.spl,tN).*cos(ppval(emb.th.spl,tN)), ...
               ppval(emb.alph.spl,tN).*sin(ppval(emb.th.spl,tN))];
    tmp_pts = (tmp_pts - repmat(b',length(tN),1))*Minv';
    % compute discrete length over original scales
    L = cumsum([0; sqrt(sum(( tmp_pts(2:end,:) - tmp_pts(1:end-1,:) ).^2,2))]);

    % refine domain using length-based importance sampling
    tN = pchip(1/max(L)*L,tN,linspace(0,1,N)');
    
elseif strcmp(smpl,'nom')
    % evaluate shape at original measure
    tmp_pts = (P - repmat(b',emb.nom.num,1))*Minv';
    L = cumsum([0; sqrt(sum(( tmp_pts(2:end,:) - tmp_pts(1:end-1,:) ).^2,2))]);
    
    % refine domain using nominal length-based importance sampling
    tN = pchip(1/max(L)*L,emb.TF.t,linspace(0,1,N)');
end

%% Evaluate shape characteristics
% save total number of points
emb.num = N;

% compute continuous curvature at updated points 
emb.curv = cont_curv(tN,emb.th.spl,emb.alph.spl,emb.th.scl,emb.alph.scl,emb.TF.Minv);
emb.nom.curv = cont_curv(emb.nom.t,emb.th.spl,emb.alph.spl,emb.th.scl,emb.alph.scl,emb.TF.Minv);

% compute continuous normal vector approximation
emb.nml     = emb.alph.scl*ppval(fnder(emb.alph.spl,1),tN).*([sin(ppval(emb.th.spl,tN)), -cos(ppval(emb.th.spl,tN))]) +...
              emb.th.scl*ppval(fnder(emb.th.spl,1),tN).*[ppval(emb.alph.spl,tN).*cos(ppval(emb.th.spl,tN)), ...
                                                         ppval(emb.alph.spl,tN).*sin(ppval(emb.th.spl,tN))];
emb.TF.nml = emb.alph.scl*ppval(fnder(emb.alph.spl,1),emb.nom.t).*([sin(ppval(emb.th.spl,emb.nom.t)), -cos(ppval(emb.th.spl,emb.nom.t))]) +...
              emb.th.scl*ppval(fnder(emb.th.spl,1),emb.nom.t).*[ppval(emb.alph.spl,emb.nom.t).*cos(ppval(emb.th.spl,emb.nom.t)), ...
                                                         ppval(emb.alph.spl,emb.nom.t).*sin(ppval(emb.th.spl,emb.nom.t))];

% evaluate shape at new N-landmarks
emb.pts = [ppval(emb.alph.spl,tN).*cos(ppval(emb.th.spl,tN)), ...
           ppval(emb.alph.spl,tN).*sin(ppval(emb.th.spl,tN))];

% save evaluated length scales
emb.t = tN;
% save evaluations of inner products
emb.nom.alph  = ppval(emb.alph.spl,emb.nom.t);
% save evaluations of angle
emb.nom.th    = ppval(emb.th.spl,emb.nom.t);

%% Transform back to original coordinates
% evaluate shape at nominal landmarks
emb.nom.pts =  ([ppval(emb.alph.spl,emb.nom.t).*cos(ppval(emb.th.spl,emb.nom.t)),...
                 ppval(emb.alph.spl,emb.nom.t).*sin(ppval(emb.th.spl,emb.nom.t))]...
                 - repmat(emb.TF.b',length(emb.nom.t),1))*emb.TF.Minv';

% transform scales back to original scale
emb.pts = (emb.pts - repmat(emb.TF.b',N,1))*emb.TF.Minv';
% compute discrete length from tesselation
emb.L = cumsum([0; sqrt(sum(( emb.pts(2:end,:) - emb.pts(1:end-1,:) ).^2,2))]);
emb.nom.L = cumsum([0; sqrt(sum(( emb.nom.pts(2:end,:) - emb.nom.pts(1:end-1,:) ).^2,2))]);
% compute angles
% emb.th.eval = unwrap(atan2(emb.pts(:,2),emb.pts(:,1))); 
emb.th.eval = ppval(emb.th.spl,tN);
% compute inner product
% emb.alph.eval = emb.pts(:,1).*cos(emb.th.eval) + emb.pts(:,2).*sin(emb.th.eval);
emb.alph.eval = ppval(emb.alph.spl,tN);

% compute unit normals
emb.nml = emb.nml*[0 -1; 1 0]*Minv'*[0 1; -1 0];
emb.nml = emb.nml./sqrt(emb.nml(:,1).^2 + emb.nml(:,2).^2);
emb.nom.nml = emb.TF.nml*[0 -1; 1 0]*Minv'*[0 1; -1 0];
emb.nom.nml = emb.nom.nml./sqrt(emb.nom.nml(:,1).^2 + emb.nom.nml(:,2).^2);