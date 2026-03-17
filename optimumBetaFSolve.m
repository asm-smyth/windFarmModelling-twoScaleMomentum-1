% --------------------------------------------------------
% Solve iteratively for optimum beta and CT by finding the roots of
% d(Cpg)/dCT = 0 (see Nishino and Smyth 2026), using fsolve
% 16-03-2026
% --------------------------------------------------------

function [beta_opt,CT_opt,CPg_opt] = optimumBetaFSolve(h0,lambda,L,Cf0,CTr,CPr,k,Cc)
chiT = 1; %placeholder 
K1 = Cf0/lambda/chiT;
K2 = h0/L/Cf0;

% Variables to solve for using fsolve:
% x = [CT0 beta0 diffBeta Psi Phi chi dChiT dChiP dSig]; 
% initial guesses:
x = [0.9 0.9 0 0 0 1 0 0 0];


fun = @(x) paramfun(x,K1,K2,lambda,Cc,k,CPr,CTr);
x = fsolve(fun,x);

beta_opt = x(2);
CT_opt = x(1);

% Define variables (from Nishino and Smyth 2026)
CPrADT = 0.5*CTr*(1+sqrt(1-CTr));
CPA = 0.5*x(1)*(1+sqrt(1-x(1)));
gama4 = (x(1)/CPA - 1)/(CTr/CPrADT - 1);
sig = sqrt( gama4 );
Eta = sig*(CPr/CPrADT -1) + 1;
CPg_opt = (x(2)^3)*(x(6)^3)*Eta*CPA;


function F = paramfun(x,K1,K2,lambda,Cc,k,CPr,CTr)
    F = [ -(3*x(5)/x(2))*(x(1)^3) + ((3*x(5)/x(2))*(1+sqrt(1-x(1))) - (3*x(3)/x(2) + x(4)))*(x(1)^2) + ((3*x(3)/x(2) + x(4))*(1+sqrt(1-x(1))) - 3/2)*x(1) + 1 + sqrt(1-x(1));
    ((x(6)^2)/K1)*(x(2)^4)/(K2*(x(2)^2) - 3*(1+K2)) - x(3);
    1 - Cc*(1-sqrt(1-x(1)))/((1 + 2*k*sqrt(pi/4/lambda))^2) - x(6);
    (x(1)*(x(6)^2)/K1 + 1)*(x(2)^3) + K2* (x(2)^2) + (-K2-1);
    x(9)*(CPr/( 0.5*CTr*(1+sqrt(1-CTr)) ) - 1)/( sqrt( (x(1)/(0.5*x(1)*(1+sqrt(1-x(1)))) - 1)/(CTr/( 0.5*CTr*(1+sqrt(1-CTr)) ) - 1) )*(CPr/( 0.5*CTr*(1+sqrt(1-CTr)) ) -1) + 1 ) + x(8)/(x(6)^3) - x(4);
    (x(3)/(x(6)^2))*x(7) - x(5);
    -Cc*( 1 - Cc*(1 - sqrt(1-x(1)))/((1 + 2*k*sqrt(pi/4/lambda))^2) )/sqrt(1-x(1))/((1 + 2*k*sqrt(pi/4/lambda))^2) - x(7);
    -Cc*(3/2)*((1 - Cc*(1 - sqrt(1-x(1)))/((1 + 2*k*sqrt(pi/4/lambda))^2) )^2)/sqrt(1-x(1))/((1 + 2*k*sqrt(pi/4/lambda))^2) - x(8);
    0.5/sqrt( (x(1)/(0.5*x(1)*(1+sqrt(1-x(1)))) - 1)/(CTr/( 0.5*CTr*(1+sqrt(1-CTr)) ) - 1) )/(CTr/(0.5*CTr*(1+sqrt(1-CTr))) - 1)/(((1 + sqrt(1-x(1)))^2)*sqrt(1-x(1))) - x(9)
    ];
end




end