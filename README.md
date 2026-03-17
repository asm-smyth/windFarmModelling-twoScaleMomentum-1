# windFarmModelling-twoScaleMomentum-1
Matlab code for calculating farm-scale induction (that is the farm-average wind speed reduction) using two-scale momentum theory and find the optimum induction given atmospheric and farm properties. Methodology based on Nishino and Smyth (2026), preprint available at: https://doi.org/10.48550/arXiv.2508.02727 

The file "main_optimumBeta" can be run directly to obtain the optimal farm performance for the input parameters specified in the file. The function "betaSolve" can be run directly to obtain the (not necessarily optimal) farm induction that results from operating at a certain  farm-average CT. 

There are two options for solving for the optimum farm induction in the "main_optimumBeta" file. One uses the matlab function fsolve, which is more robust but requires the Matlab Optimization Toolbox. The other uses simple manual iteration. Either option can be commented/uncommented in the script. 

The inputs that must be specified are:
L - the farm length in the wind direction [m].
Cf0 - the ground friction coefficient.
h0 - the atmospheric boundary layer height.
CTr - the rated thrust coefficient of the turbines in the farm.
CPr - the rated power coefficient of the turbines in the farm.
lambda - the density of turbines in the farm, which for a homogeneous farm is given by lambda = pi/(4(s/D)^2) where D is the turbine diameter and s is the distance between consecutive turbine rows.

Furthermore, two empirical constants are required (both defined in Nishino and Smyth 2026):
k - wake expansion factor.
Cc - layout coefficient.

The key outputs are:
beta_opt - the farm speed reduction factor beta which corresponds to maximum farm power generation.
CT_opt - the farm-average CT that results in the optimum farm speed reduction factor.
CPg_opt - the farm efficiency that results from operating at optimal induction.

