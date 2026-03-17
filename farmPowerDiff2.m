% --------------------------------------------------------
% Function to solve cubic equation for CT (equation 29 in Nishino and 
% Smyth, 2026), as a step in iteratively solving for optimum beta and the 
% CT that results in optimum beta. The function also calculates the global 
% farm efficiency CPg for the obtained values of CT and beta.
% 16-03-2026
% --------------------------------------------------------

function [CT,CPg,beta] = farmPowerDiff2(beta0,diffBeta,CT0,lambda,h0,L,Cf0,CTr,CPr,k,Cc,Phi)
chi = 1 - Cc*(1-sqrt(1-CT0))/((1 + 2*k*sqrt(pi/4/lambda))^2);
chiT = chi^2;
chiP = chi^3;

        % Accounting for non-ideal turbine
        CPrA = 0.5*CTr*(1+sqrt(1-CTr));

        CPA = 0.5*CT0*(1+sqrt(1-CT0));
        gama4 = (CT0/CPA - 1)/(CTr/CPrA - 1);
        sig = sqrt( gama4 );
        gama5 = CTr/CPrA - 1;
        gama6 = ((1 + sqrt(1-CT0))^2)*sqrt(1-CT0);
        dSig = 0.5/sig/gama5/gama6;
        Eta = sig*(CPr/CPrA -1) + 1;
        dEta = dSig*(CPr/CPrA - 1);

        % Accounting for non-ideal farm layout
        dChiP = -Cc*(3/2)*((1 - Cc*(1 - sqrt(1-CT0))/((1 + 2*k*sqrt(pi/4/lambda))^2) )^2)/sqrt(1-CT0)/((1 + 2*k*sqrt(pi/4/lambda))^2);
       
        % Assemble coefficients for cubic equation
        Psi = dEta/Eta + dChiP/chiP;
        gama1 = 3*Phi/beta0;
        gama2 = 3*diffBeta/beta0 + Psi; 
        gama3 = sqrt(1-CT0);

    % Cubic solution of Eq. 29 from Nishino and Smyth (2026)
    a = -gama1;
    b = gama1*(1+gama3) - gama2;
    c = gama2*(1+gama3) - 3/2;
    d = 1 + gama3;
    p = [a b c d]; CTRoots = roots(p);
    
    % Identify correct root (real value in range 0 to 1)
    if sum(imag(CTRoots)==0) >= 2
        CT = min((imag(CTRoots)==0).*CTRoots);
    else
        CT = max((imag(CTRoots)==0).*CTRoots);
    end

% Calculate global farm efficiency CPg 
CPg = (beta0^3)*chiP*Eta*CPA;

% Calculate beta corresponding to the value of CT found from the cubic solution
beta = betaSolve(CT,h0,lambda,L,Cf0,k,Cc);

end