% Random hypercube corner sweeps with linear inequality
% WARNING: The nominal design must be a feasible non-corner point
close all; clearvars; rng(48);

%% Algorithm parameters (things to modify)
% number of sweeps
N = 10;
% number of points over sweep
T = 20; t = linspace(0,1,T);
% dimension of domain
m = 3;

%% Convex domain definition
% nominal
x0 = zeros(m,1);
% some constraint choices
% A = ones(1,m); b = 1; % one active constraint 
A = [ones(1,m);-ones(1,m)]; b = [1; 1]; % two active constraints
% A = ones(1,m); b = 10; % non-active constraint
% some boundary choices
% xu = linspace(1,10,m)'; xl = -0.5*ones(m,1);
xu = ones(m,1); xl = -ones(m,1);

%% Hypercube transformations
M    = diag(2./(xu-xl)); a    = -(M*xl + ones(m,1)); % i.e., x  = M*x0   + a    in [-1,1]_m
Minv = diag((xu-xl)./2); ainv = -Minv*a;             % i.e., x0 = Minv*x + ainv in [xl,xu]_m

%% Determine corners for sweeps
% check for feasible baseline 
if max(A*x0 - b) <= 0 && min(x0 <= xu) == 1 && min(x0 >= xl)
    disp('Baseline design is feasible...')
else
    disp('WARNING: Baseline design is not feasible')
end

% precondition corner matrix and counter
XC = zeros(m,2*N); i = 1;
% linprog() options
options = optimset('linprog'); options.Display = 'off';
% find initial corner point
p = 2*rand(m,1) - 1; XC(:,i) = linprog(p',A*Minv,b - A*ainv,[],[],-ones(m,1),ones(m,1),options);
% use the oposite direction to find another corner
i = 2; XC(:,i) = linprog(-p',A*Minv,b - A*ainv,[],[],-ones(m,1),ones(m,1),options);
% loop to find additional corners if necessary
if N > 1 && nchoosek(i,2) < N
while nchoosek(i,2) < N && i < 2^m - (m+1)
    % pick new random direction
    p = 2*rand(m,1) - 1; [c, fval] = linprog(p',A*Minv,b - A*ainv,[],[],-ones(m,1),ones(m,1),options);
    % criteria for determining if point is a new corner (only valid for non-zero corners, good thing we scaled!)
    if min( sum((XC(:,1:i) - repmat(c,1,i)).^2,1) ) ~= 0
        i = i + 1;
        disp(['New corner found... ',num2str(nchoosek(i,2)/N*100),'% complete (',num2str(nchoosek(i,2)),'/',num2str(N),' possible sweeps)'])
        XC(:,i) = c;
    end
end
else
    disp(['Complete with ',num2str(nchoosek(i,2)),' possible sweep(s) of requested ',num2str(N),'.'])
end

%% Build sweeps from found corners
% resize to found corners (remove any remaining preconditioning)
XC = XC(:,1:i);
% compute corner sweeps from combinations
cind = nchoosek(linspace(1,i,i),2);
% sort on largest sweeps
[~,sind] = sortrows(sum((XC(:,cind(:,1))' - XC(:,cind(:,2))').^2,2),'descend');
% compute the first nchoosek(i,2) or N largest sweeps from the corner combinations
if nchoosek(i,2) <= N
    S  = kron(XC(:,cind(sind(1:nchoosek(i,2)),1)),(1-t))' + kron(XC(:,cind(sind(1:nchoosek(i,2)),2)),t)';
    % scale back to original units
    S0 = S*Minv + repmat(ainv',nchoosek(i,2)*T,1);
elseif nchoosek(i,2) > N
    S  = kron(XC(:,cind(sind(1:N),1)),(1-t))' + kron(XC(:,cind(sind(1:N),2)),t)';
    % scale back to original units
    S0 = S*Minv + repmat(ainv',N*T,1);
end

%% Build remaining sweeps, if necessary, from random directions through nominal
% determine if number of sweeps has not been met (fill in with random directions)
Nrmg = N-nchoosek(i,2);
if Nrmg > 0
    disp(['Filling in ',num2str(Nrmg),' remaining sweeps with random directions...']);
    % sample random directions
    P = (2*rand(m,Nrmg) - 1); P = P./sqrt(sum(P.^2,1));
    % scale and concatenate constraints
    AA = [A*Minv; eye(m); -eye(m)]; bb = [b - A*ainv; ones(2*m,1)]; 
    % compute extents   
    tt = (bb - AA*repmat(M*x0 + a,1, Nrmg))./(AA*P); ext = tt;
    % take first two minimum magnitude extents for each random direction
    % sort on magnitudes
    [~,tind] = sort( abs(tt) );
    % loop through columns to build the two sorted-extent vectors
    for j = 1:N-nchoosek(i,2)
        ext(:,j) = tt(tind(:,j),j);
        % take opposite side (sign must change) as second extent
        if sign(ext(1,j)) ~= 0 % non-active nominal
            extsgn = ext((sign(ext(1,j)).*sign(ext(:,j)) == -1),j);
        else                   % active nominal
            extsgn = ext((sign(ext(1,j)) - sign(ext(:,j)) ~= 0),j);
        end
        ext(2,j) = extsgn(1,1);
    end  
    % compute boundary points
    X  = [( P.*ext(1,:) + repmat(M*x0 + a,1,Nrmg) )';...
          ( P.*ext(2,:) + repmat(M*x0 + a,1,Nrmg) )'];
    % sweep over boundary points and append
    Srnd = kron(X(1:Nrmg,:)',(1-t))' + kron(X(Nrmg+1:end,:)',t)';
    S  = [S; Srnd];
    % scale to original units and append
    S0 = [S0; Srnd*Minv + repmat(ainv',Nrmg*T,1)];
    %% Check feasibility
    if max(max(AA*S' - bb)) > eps
        disp(['WARNING: Some point(s) violate constraint(s). Worst case scaled-violation = ',num2str(max(max(AA*S' - bb)))]);
    end
end
%% Visualize
% 2-dimensional visualization
if m == 2
    close all; scatter(XC(1,:),XC(2,:),'filled','k'); hold on;
    scatter(S(:,1),S(:,2),'k.');
    scatter(X(:,1),X(:,2),'k');    
end
% 3-dimensional visualization
if m == 3
    close all; scatter3(XC(1,:),XC(2,:),XC(3,:),'filled','k'); hold on;
    scatter3(S(:,1),S(:,2),S(:,3),'k.');
    scatter3(X(:,1),X(:,2),X(:,3),'k');    
end
    
%% Write results
csvwrite('sweeps.csv',S0)