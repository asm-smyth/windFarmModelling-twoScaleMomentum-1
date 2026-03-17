
% --------------------------------------------------------
% Solve iteratively for optimum beta and CT by finding the roots of
% d(Cpg)/dCT = 0 (see Nishino and Smyth 2026)
% 16-03-2026
% --------------------------------------------------------

function [beta_opt,CT_opt,CPg_opt] = optimumBeta(h0,lambda,L,Cf0,CTr,CPr,k,Cc)

% initial guesses
CT0 = 0.9;
beta0 = 0.9;

% Error functions
err1 = 1;
err2 = 1;
maxErr = 0.0001;

% Iterate until convergence
ii = 0;
while err1 > maxErr || err2 > maxErr
    ii=ii+1;

    % Find beta and d(beta)/dCT from momentum equation
    [diffBeta,beta,Phi] = twoScaleMom2(beta0,CT0,lambda,h0,L,Cf0,k,Cc); 
    % Find CT by solving d(CPg)/dCT = 0
    [CT,CPg,beta] = farmPowerDiff2(beta0,diffBeta,CT0,lambda,h0,L,Cf0,CTr,CPr,k,Cc,Phi); 

    % Evaluate step error
    err1 = abs(CT-CT0)/CT0;
    err2 = abs(beta-beta0)/beta0;
    CT0 = CT;
    beta0 = beta;
end

% Assign outputs
CT_opt = CT;
beta_opt = beta;
CPg_opt = CPg;


end