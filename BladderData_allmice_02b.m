%% Bladder Cancer Data for all mice
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Ultrasound data over time for each of the 55 experimental 
% mice from 4 experimental cohorts: untreated, Gem, OT-1, Gem+OT-1.


%% Volume data in terms of mm^3 

% All ultrasound data
Timedata = [6; 9; 13; 16; 20];
ultrasound_alldata = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_ODE_VC_GitHub/ultrasound_data.xlsx','Sheet','data','Range','B2:BI7'))';
% untreated, gem mono, ot1 mono, gem+ot1 combo--each have 15 columns of 5 time points

% when tumor volume hits 250, that is not an actual data point--they sacked the mouse
for i = 1:length(ultrasound_alldata(:,1))
    index = find(ultrasound_alldata(i,:)==250);
    ultrasound_alldata(i,index) = NaN(length(index),1);
end

%% Sort each mouse into correct cohort

% Lowerbound and upperbound of indices for when cohorts start and end
lb = [1,16,31,46]; 
ub = [14,29,42,60];

num_cohorts = length(lb);

exp_individual = nan(num_cohorts,length(Timedata),max(ub-lb)+1); % cohorts x times x number of mice
for cohort = 1:num_cohorts
    exp_individual(cohort,:,1:ub(cohort)-lb(cohort)+1) = ultrasound_alldata(lb(cohort):ub(cohort),1:length(Timedata))';
end

%% Get times data for each mouse

times_all = repmat(reshape(Timedata,1,[],1), num_cohorts, 1, max(ub-lb)+1); % cohorts x times x number of mice
times_all(isnan(exp_individual)) = NaN;


