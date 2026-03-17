% 13-03-2026
% Author: Dr Amanda Smyth (amanda.smyth@eng.ox.ac.uk)
%
% This Matlab script calculates the optimal wind farm induction using the 
% method outlined in the paper by Nishino and Smyth (2026, preprint: 
% https://doi.org/10.48550/arXiv.2508.02727), for given farm and 
% atmospheric parameters. The optimum induction condition is the 
% farm-average wind speed reduction for which the farm produces the most 
% power. 
% 
% The script gives the optimum farm-average induction in terms of the wind 
% speedreduction paramter beta = U_F/U_F0, where U_F and U_F0 are the 
% average wind speeds in the farm area with and without the farm present 
% respectively (that is, beta represents the reduction in area-average wind 
% speed due to the presence of the farm).The script also gives the 
% farm-average turbine thrust coefficient (CT) that results in the optimal 
% beta.
%
% The beta that results from operating at any given CT can also be obtained
% by simply calling the function "betaSolve" for the desired input 
% parameters.
%
%--------------------------------------------------------------------------

clear;

% ---------------------------- **MODEL INPUTS** ---------------------------
% Specify farm and atmospheric properties ---------------------------------

L = 20000; % Farm length in the wind direction [m]
Cf0 = 0.002; % Ground friction coefficient
h0 = (L*Cf0)*linspace(10,100,10); % Boundary layer height [m].
% NOTE: either a single value or an array can be specified for h0.

% Specify turbine properties ----------------------------------------------

CTr = 0.8; % Rated turbine thrust coefficient
CPr = 0.489; % Rated turbine power coefficient
lambda = Cf0*linspace(10,20,20); % Farm density
% NOTE: either a single value or an array can be specified for lambda.

% Specify empricial constants for the analytical models -------------------
% See Nishino and Smyth (2026) for meaning and discussion of each

k = 0.05; % Wake expansion factor
Cc = 0.14; % Layout coefficient 

% ---------------------------**END OF MODEL INPUTS**-----------------------
% -------------------------------------------------------------------------

nb = length(lambda);
mb = length(h0);

% -------------------------------------------------------------------------
% ------------------------------**OUTPUTS:**-------------------------------
% beta_opt -- [mb,nb] array with optimal beta values
% CT_opt -- [mb,nb] array with farm-average thrust coefficient CT that results in optimum beta
% CPg_opt -- [mb,nb] array with global farm efficiency CPg that results from operating at optimum beta
% ---------------------------**END OF MODEL OUTPUTS**----------------------
% -------------------------------------------------------------------------

% Calculate optimum beta iteratively --------------------------------------

beta_opt = zeros(mb,nb); CT_opt = beta_opt; CPg_opt = beta_opt;
for j=1:nb
    for i=1:mb
% ------------TWO FUNCTION OPTIONS, UNCOMMENT DESIRED OPTION---------------
% ** OPTION 1 ** : manual iteration (less robust)
        % [beta_opt(i,j),CT_opt(i,j),CPg_opt(i,j)] = optimumBeta(h0(i),lambda(j),L,Cf0,CTr,CPr,k,Cc);
% ** OPTION 2** : iteration using fsolve (robust, requires the Matlab Optimization Toolbox)
        [beta_opt(i,j),CT_opt(i,j),CPg_opt(i,j)] = optimumBetaFSolve(h0(i),lambda(j),L,Cf0,CTr,CPr,k,Cc);
% -------------------------------------------------------------------------
    end
end

disp('Iterative solution complete')


% -------------------------------------------------------------------------
% ------------------------------**PLOT RESULTS:**--------------------------
% As an example, the figures below plot one ABL height result against 
% effective array density, and one array density result against effective 
% ABL height.
% -------------------------------------------------------------------------

fn=10;

figure(fn); hold on
plot(lambda/Cf0,beta_opt(round(mb/2),:),'r',lambda/Cf0,CT_opt(round(mb/2),:),'b',lambda/Cf0,CPg_opt(round(mb/2),:),'k')
xlabel('\lambda/C_f_0')
title(strcat('h_0/LC_f_0 = ',num2str(h0(round(mb/2))/L/Cf0)))
legend('\beta_{(opt)}','C_{T(opt)}','C_{PG(max)}')
grid on
hold off; fn = fn+1;

figure(fn); hold on
plot(h0/Cf0/L,beta_opt(:,round(nb/2)),'r',h0/Cf0/L,CT_opt(:,round(nb/2)),'b',h0/Cf0/L,CPg_opt(:,round(nb/2)),'k')
xlabel('h_0/LC_f_0')
title(strcat('\lambda/C_f_0 = ',num2str(lambda(round(nb/2))/Cf0)))
legend('\beta_{(opt)}','C_{T(opt)}','C_{PG(max)}')
grid on
hold off; fn = fn+1;

% -------------------------------------------------------------------------
% END OF SCRIPT