% --------------------------------------------------------
% Solving for the farm-average wind speed reduction (beta) that results
% from operating at a given farm-average thrust coefficient (CT), by 
% solving the cubic Eq. 21 in Nishino and Smyth (2026).
% 17-03-2026
% --------------------------------------------------------

function [beta] = betaSolve(CT,h0,lambda,L,Cf0,k,Cc)

chi = 1 - Cc*(1-sqrt(1-CT))/((1 + 2*k*sqrt(pi/4/lambda))^2);
chiT = chi^2;

K1 = Cf0/lambda/chiT;
K2 = h0/L/Cf0;
betaRoots = roots([(CT/K1 + 1) K2 0 (-K2-1)]);
beta = max((imag(betaRoots)==0).*betaRoots);

end