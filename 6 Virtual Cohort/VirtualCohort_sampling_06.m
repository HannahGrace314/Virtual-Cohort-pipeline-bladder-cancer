%% Virtual cohort sampling - accept or reject method
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Choose which parameters to vary for LHS sampling using
% chooseparam and then code generates n sample subjects as well as n 
% initial conditions based on the distribution of data at day 6.
% Virtual subjects that produce simulations within 1, 2, ..., nstd 
% standard deviations from the mean ultrasound data for each time point 
% and for each treatment cohort (from 1,..., num_samplecohort) are marked
% with a 1 in accepted_withinstd_all, which has dimensions of n x nstd.

% Uses: BladderData_shortsummary_02.m, BladderParam_03.m, BladderFunc_01.m,
% and sample_param_03c.m
 

%% Customize sampling

% Choose Parameters ((0) - parameter fixed, (1) - parameter varies): 
chooseparam = [0;       % 1  - pc
               0;       % 2  - cmax
               0;       % 3  - ktc
               0;       % 4  - kgc
               1;       % 5  - nct
               0;       % 6  - smt
               0;       % 7  - kgt
               0;       % 8  - dt
               0;       % 9  - T0
               1;       % 10 - rcm
               1;       % 11 - kgm
               0;       % 12 - stm
               0;       % 13 - dm
               0;       % 14 - M0
               0;       % 15 - dg
               0];      % 16 - km

% Number of sample parameter sets:
    n = 500000; 
                          
% Test acceptance up to 1, 2, ..., nstd from mean:
    nstd = 4; 

% Number of cohorts to use in acceptance (advised to leave one out for validation):
    num_samplecohort = 3; % ie first 3 cohorts


%% Sampling parameter sets

% Specify range for parameters (also specifies BCFunc):
BladderParam_03

% Sample sets using Latin Hypercube Sampling:
sample_param_03c


%% Determine various initial conditions using distribution of data at day 6 
BladderData_shortsummary_02

% Normal distribution for total tumor volume at day 6
pd = makedist('Normal','mu',ultrasound_meandata(1,1),'sigma',ultrasound_stddata(1,1));
tpd = truncate(pd,max(0, ultrasound_meandata(1,1) - nstd*ultrasound_stddata(1,1)),nstd*ultrasound_stddata(1,1));
ic_ultrasound_n = random(tpd,n,1);

% Normally distributed initial condition at day 6
initialcondition_n = ic_ultrasound_n*[initialpercent,0];


%% Calculate numerical simulations, error, and if sets within STD threshold 
timerange = (Timedata(1):1:23)'; % Evaluate_NumSim_03g.m uses a time step of 1
collect_numsim = nan(num_samplecohort,4,length(timerange), n);
accepted_withinstd = zeros(n,nstd,num_samplecohort);

for cohort = 1:num_samplecohort
    % Treatment Modes:
    % 1 - untreated
    % 2 - Gem
    % 3 - OT-1 
    % 4 - Gem+OT-1 

    % Load data
    dose = cohort_dose(cohort,:);
    ultrasound_mean = ultrasound_meandata(:,cohort);
    ultrasound_std = ultrasound_stddata(:,cohort);
    
    %----------------------- Numerical simulations ------------------------
    nd = length(Timedata);
    totalvolume_numsim = zeros(nd,n);
    parfor mouse = 1:n
        mouse
        pull_numsim = Evaluate_NumSim_03g(BCFunc,Psample(mouse,:),[Timedata(1),23],initialcondition_n(mouse,:),dose);
        X = pull_numsim{1};
        T = pull_numsim{2};
        numsim(:,:,mouse) = [X, sum(X,2)]';
        totalvolume_numsim(:,mouse) = sum(X(ismember(T, Timedata),:),2); % total tumor at data collection points
    end

    %---------------- Keep simulations for future plotting ----------------
    collect_numsim(cohort,:,:,:) = numsim;
    
    %------------------------- Calculate error ----------------------------
    voldata_matrix = ones(nd,n).*ultrasound_mean;
    stddata_matrix = ones(nd,n).*ultrasound_std;
    
    error_ateachday = abs(totalvolume_numsim - voldata_matrix);
    
    % ------- Keep parameter sets within a specified STD threshold --------
    for i = 1:nstd
        withinstd_ateachday = error_ateachday <= i * stddata_matrix;
        accepted_withinstd(:,i,cohort) = prod(withinstd_ateachday);
    end

end


%% Accept/reject parameter sets

accepted_withinstd_all = logical(prod(accepted_withinstd,3)); 

