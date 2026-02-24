%% Least Squares Error Cost Function for model compared to ultrasound data
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Calculates the least squares error of the total tumor volume 
% computed by the model versus the ultrasound data. Used in computing 
% profile likelihoods (Figure 8) in identifiability_main_04.m.

function error = ModelCost_ultrasound_04b(BCFunc,params, times,data,timerange,initialcondition,dose)

% Evaluate ODE:
numsim = Evaluate_NumSim_03g(BCFunc,params,timerange,initialcondition,dose);
T = numsim{2};
X = sum(numsim{1},2);

% Pull out data times:
x = X(ismember(T, times), :);
t = nonzeros(ismember(T, times).*T);

% Least Squares Error:
error = sum((data - x).^2); 


end
