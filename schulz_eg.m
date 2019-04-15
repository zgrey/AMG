% Riemannian view on shape optimization [2014]
% V. Schulz numerical example
close all; clear variables; clc;
addpath './Euclidean_tools/Euclidean_Shapes/';

%% Things to modify
N      = 101;   % Number of discrete shape landmarks
mu     = 2;     % function parameter (see [Schulz, 2014])
method = 'NM';  % optimization method (NM or SD)
maxi   = 20;    % max number of iterations

% initial shape
theta = linspace(0,2*pi,N)';
% Schulz's initial shape
c0 = 1/2*[cos(theta)-0.15*abs(1-sin(2*theta)).*cos(theta),sin(theta)-0.15*abs(1-cos(2*theta)).*cos(theta)];
% circle
% c0 = [cos(theta), sin(theta)];
% tall ellipse
% a = 1; b = 2; c0 =[a*cos(theta),b*sin(theta)];
% Airfoil
% load turbine_airfoil.mat; x0 = [mean(PP(:,1)),mean(PP(:,2))]; c0 = PP - x0;

% plot nominal shape
[c0,n0] = embrep3(c0,size(c0,1),'uni');
h = plot(c0(:,1),c0(:,2),'o-'); axis equal; hold on; title 'Initial Shape'
% throw out last repeated point
c0 = c0(1:end-1,:); n0 = n0(1:end-1,:);
quiver(c0(:,1),c0(:,2),n0(:,1),n0(:,2),'color',h.Color);

% functional approximation
f = @(c) 1/(2*mu)*sum( ( unwrap(atan2(mu*c(2:end,2),c(2:end,1)))-unwrap(atan2(mu*c(1:end-1,2),c(1:end-1,1))) )...
    .*( (c(2:end,1).^2 + mu^2*c(2:end,2).^2).^2/4 + (c(1:end-1,1).^2 + mu^2*c(1:end-1,2).^2).^2/4 ...
    - (c(2:end,1).^2 + mu^2*c(2:end,2).^2)/2 - (c(1:end-1,1).^2 + mu^2*c(1:end-1,2).^2)/2) );

% shape distance approximation
d = @(c) 1/(2*mu)*sum( ( unwrap(atan2(mu*c(2:end,2),c(2:end,1)))-unwrap(atan2(mu*c(1:end-1,2),c(1:end-1,1))) )...
    .*(abs(sqrt(c(2:end,1).^2 + mu^2*c(2:end,2).^2)-1) + ...
    abs(sqrt(c(1:end-1,1).^2 + mu^2*c(1:end-1,2).^2)-1)) );

%% Schulz shape optimizatoin example
figure; h = plot([c0(:,1);c0(1,1)],[c0(:,2);c0(1,2)]); axis equal; hold on; title Iterates;
% quiver(c0(:,1),c0(:,2),n0(:,1),n0(:,2),'color',h.Color);
% preallocate for speed
a_opt = zeros(maxi,1); CI = zeros(maxi,1); F = zeros(maxi,1); D = zeros(maxi,1);
% initialize solvers and counter
c = c0; n = n0; gradf = ones(N,1); i = 1;
while norm(gradf) > 1e-5 && i <= maxi && d(c) > 1e-8
    % functional
    F(i) = f(c);
    % distance metric
    D(i) = d(c);
    
    % shape gradient per functional defined on pp. 496 equation (4.1)
    gradf = c(:,1).^2 + mu^2*c(:,2).^2 - 1;
    % compute Hessian at discrete landmarks per paper pp. 496 equation (4.3)
    s = c./sqrt(c(:,1).^2 + mu^2*c(:,2).^2);
    Hessf = 2*sqrt(s(:,1).^2 + mu^4*s(:,2).^2);
    
    if strcmp(method,'NM')
        % compute scalar-valued function for retraction map per Newton Method
        R = -gradf./Hessf;
    elseif strcmp(method,'SD')
        % compute scalar-valued function for retraction map per steepest-descent
        R = -gradf;
        
    end
    
    % exact line search
    a_opt(i) = fminbnd(@(alpha) f(c + alpha*R.*n),0.1,1);   
    
    % update shape via the retraction map pp. 495 
    c = c + a_opt(i)*R.*n;
    
    % recompute shape with uniform measure
    [c,n] = embrep3(c,size(c,1)+1,'curv'); 
    % throw out last repeated point
    c = c(1:end-1,:); n = n(1:end-1,:);
    
    % convergence criteria
    CI(i) = norm(gradf);
%     CI(i) = d(c);
    
    h = plot([c(:,1);c(1,1)],[c(:,2);c(1,2)]); hold on;
    quiver(c(:,1),c(:,2),n(:,1),n(:,2),'color',h.Color);

    i = i + 1;
end
i = i - 1;
% convergence plots
figure; plot([c(:,1);c(1,1)],[c(:,2);c(1,2)],'k-o','linewidth',1); hold on;
quiver(c(:,1),c(:,2),n(:,1),n(:,2),1,'k');
axis equal; title 'Final Shape';
figure; semilogy(1:i,CI(1:i),'k-o','linewidth',2); grid on; title Convergence; hold on; 
ylabel '||grad(f)||'; xlabel 'number of iterations'; grid on;

% results from paper
semilogy([1:5],[0.9222,0.1382,0.3571e-2,0.8187e-5,0.1736e-9],'b-o','linewidth',2);
dSD = [0.9222e00; 0.2137e00; 0.6174e-01;0.1730e-01;0.4888e-02;0.1448e-02 ;0.4404e-03 ;0.1300e-03 ;0.4065e-04 ;0.1228e-04 ;0.3909e-05 ;0.1198e-05 ;0.4026e-06 ;0.1171e-06 ;0.3650e-07 ;0.8370e-08 ;0.3041e-08];
semilogy([1:17],dSD,'r-o','linewidth',2);
legend('This implementation','d(c) from Schulz Paper (NM)','d(c) from Schulz Paper (SD)');