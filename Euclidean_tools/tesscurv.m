% Discrete curvature based on tangent orientation

function [kappa] = tesscurv(CC)

tbwd = [CC(end,:) - CC(end-1,:) ; CC(2:end-1,:) - CC(1:end-2,:)];
tau  = (CC(2:end,:) - CC(1:end-1,:));

% tangent vector angles
ang = acos(sum(tau.*tbwd,2) ./ ( sqrt(sum(tau.*tau,2)).*sqrt(sum(tbwd.*tbwd,2)) ));
% discrete curvature
kappa = 2*ang ./ (sqrt(sum(tau.*tau,2)) + sqrt(sum(tbwd.*tbwd,2))); kappa = [kappa;kappa(1)];