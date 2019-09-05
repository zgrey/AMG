% Non-uniform FFT

function [ff,xx] = NUFFT(x,f,NN)

N = length(f);
xx = linspace (-1, 1, NN); 

% construct DFT matrix 
k = [0:floor(N/2), -floor(N/2)+1:-1];
A = exp(1i * x * k * pi / 2); 

% solve y = a * f and obtain f, i.e. nuDFT of y. 
% exactly, to stabilize the solution, solving L2 penalized 
% simeq: a' * y = (a' * a + alpha^2 * I) * f, 
% and its solution minimize (y - a * f)^2 + alpha^2 * |f|^2 
alpha = 2; 
z = (A' * A + alpha^2 * eye(size(A,2))); 
c = z \ (A' * f); 

% iDFT 
a0 = exp(1i * xx' * k * pi / 2); ff = a0 * c;