% --------------------------------------------------------
% Function to calculate the gradient of farm-average wind reduction (beta)
% with respect to CT based on two-scale momentum theory, as a step in 
% iteratively solving for optimum beta and the CT that results in optimum 
% beta. This function solvies Eq. 23 in Nishino and Smyth (2026). It also
% calculates the non-ideal turbine layout parameter Phi (Eq. 24).
% 17-03-2026
% --------------------------------------------------------


function [diffBeta,beta,Phi] = twoScaleMom2(beta0,CT0,lambda,h0,L,Cf0,k,Cc)

chi = 1 - Cc*(1-sqrt(1-CT0))/((1 + 2*k*sqrt(pi/4/lambda))^2);
chiT = chi^2;

    beta = betaSolve(CT0,h0,lambda,L,Cf0,k,Cc);

    K1 = Cf0/lambda/chiT;
    K2 = h0/L/Cf0;
    diffBeta = (1/K1)*(beta0^4)/(K2*(beta0^2) - 3*(1+K2));
    
    % Accounting for non-ideal layout effects
    dChiT = -Cc*( 1 - Cc*(1 - sqrt(1-CT0))/((1 + 2*k*sqrt(pi/4/lambda))^2) )/sqrt(1-CT0)/((1 + 2*k*sqrt(pi/4/lambda))^2);
    Phi = (diffBeta/chiT)*dChiT;


end