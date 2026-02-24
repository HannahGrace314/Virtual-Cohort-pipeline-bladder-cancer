%% Parameter fitting 
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Fits parameters hiearchically two cohorts at a time
% depending on if the chosen_cohort is gem-treated (gem and gem+ot1 cohorts)
% or not gem-treated (untreated and ot1 cohorts). First, code finds an 
% initial guess for gradient descent by using latin hypercube sampling 
% and then choosing the set of the lowest error. Afterwards, gradient descent 
% is performed to find the best fit based on data from those chosen two cohorts.

% Uses: BladderData_shortsummary_02.m, BladderParam_03.m, BladderFunc_01.m,
% ChooseParamFunc_03e.m, sample_param_03c.m, fminsearchbnd.m

clear 

n = 10000; % number of sample parameters from range

%% Choose Cohorts for Fitting

%chosen_cohort = [1,3]; % untreated, ot1
chosen_cohort = [2,4]; % gem, gem/ot-1

% Cohorts:
% 1 - untreated
% 2 - Gemcitabine mono
% 3 - OT-1 mono 
% 4 - Gem/OT-1 combo 

%% Choose Parameters to Fit

chooseparam = [0;       % 1  - pc
               0;       % 2  - cmax
               0;       % 3  - ktc
               1;       % 4  - kgc
               0;       % 5  - ntc
               0;       % 6  - smt
               1;       % 7  - kgt
               0;       % 8  - dt
               0;       % 9  - T0
               0;       % 10 - rcm
               1;       % 11 - kgm
               0;       % 12 - stm
               0;       % 13 - dm
               0;       % 14 - M0
               0;       % 15 - dg
               0];      % 16 - km

%% Sample parameters

% File loads BCFunc = BladderFunc_01.m, ranges, and values for parameters
BladderParam_03

% File samples parameters with Latin Hypercube Sampling
sample_param_03c

%% Load Data for input into cost function 
BladderData_shortsummary_02

num_chosen_cohorts = length(chosen_cohort);
all_data = cell(num_chosen_cohorts,1);
all_dose = cell(num_chosen_cohorts,1);
for j = 1:num_chosen_cohorts
    cohort = chosen_cohort(j);

    % Load data
    cell_voldata = squeeze(celltype_meandata(:,cohort,:));
    initialcondition = [cell_voldata(1,:), 0]; 
  
    num_datavariables = length(cell_voldata(1,:));
    times = [Timedata, [NaN; NaN; NaN; 16; 20], [NaN; NaN; NaN; 16; 20]];
    data = NaN(length(Timedata),2,num_datavariables);
    for i = 1:num_datavariables
        data(:,:,i) = [times(:,i),cell_voldata(:,i)];
    end
    
    all_data{j} = data;
    all_dose{j} = cohort_dose(cohort,:);
end


%% Calculate the relative error of each sampled set compared to data 
error_cohort = zeros(n,num_chosen_cohorts);

for j = 1:num_chosen_cohorts
    for i = 1:n
        error_cohort(i,j) = ModelCost_histologyflow_03d(BCFunc, Psample(i,:), [6,23], all_data{j}, initialcondition, all_dose{j});
    end
end


%% Sampled parameter set of minimum total relative error

error_total = sum(error_cohort,2);
minP_total = Psample(error_total == min(error_total),:);


%% Initial guess for fixed and varying parameters

% Starting parameter values:
fixed = minP_total'; 

% Parameters we are varying:
indexchooseparam = find(chooseparam == 1);   % index
initialguess = nonzeros(chooseparam.*fixed); % inital guess
paramrange = mpr(indexchooseparam,:);        % range


%% Parameter estimation 
% Gradient descent - estimate the model parameters:
%   if exitflag = 1,  the function converged to a solution
%   if exitflag = 0,  the number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.MaxFunEvals.
%   if exitflag = -1, the algorithm was terminated by the output function.

[paramests, fval, exitflag] = fminsearchbnd(@(p) ...
    allcohort_error(BCFunc, ChooseParamFunc_03e(fixed,p,indexchooseparam), [6,23], all_data, initialcondition, all_dose, num_chosen_cohorts),...
    initialguess,paramrange(:,1),paramrange(:,2),optimset('Display','iter','MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-7,'TolFun',1e-7));

paramests = abs(paramests);
gradientdescent_estimate = ChooseParamFunc_03e(fixed,paramests,indexchooseparam);


%% Function to calculate the sum of error across all cohorts evaluated
function allcohort_error = allcohort_error(BCFunc, params, timerange, data, initialcondition, dose, num_chosen_cohorts)

    allcohort_error = 0;
    for j = 1:num_chosen_cohorts
        allcohort_error = allcohort_error + ModelCost_histologyflow_03d(BCFunc, params, timerange, data{j}, initialcondition, dose{j});
    end


end



