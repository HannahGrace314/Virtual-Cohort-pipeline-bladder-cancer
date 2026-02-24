%% Practical Identifiability with Profile Likelihoods
% Original Source: Marisa Eisenberg (marisae@umich.edu) 
% Modified to evaluate parameters over multiple experimental cohorts for 
% Anderson et al. (2026) -- "A virtual cohort framework with applications 
% to adoptive cell therapy in bladder cancer"

% Description: Produces profile likelihoods for chosen parameters
% (chooseparam) based on ultrasound data from multiple experimental
% cohorts. Plots profile likelihoods and 95% confidence threshold as well.

% Uses: ModelCost_ultrasound_04b.m, ProfLike_04c.m, BladderParam_03.m, 
% BladderFunc_01.m, and BladderData_shortsummary_02

%% Profile likelihood specifications

% Number of points at which to sample profile likelihood:
numpoints = 200; % increasing numpoints results in a smoother profile, but will make code take longer

% Sampling range for profile likelihood:
profrange = 0.7; % profrange determines the range profile likelihoods sampled over. 
                 % profrange = 1 means that parameter ranges from 0 to 2*estimated value

% Confidence level for threshold:
alpha = 0.05; % usually to 0.05 so you can get a 95% confidence interval 


%% Choose parameters to vary
% (0) excludes or (1) includes parameters

chooseparam = [0;       % 1  - pc
               0;       % 2  - cmax
               1;       % 3  - ktc
               0;       % 4  - kgc
               0;       % 5  - ntc
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

                
%% Load Data
BladderData_shortsummary_02 % has ultrasound_meandata, Timedata, celltype_meandata (for initial condition), and cohort_dose
num_cohorts = size(ultrasound_meandata,2); % number of treatment cohorts


%% Initial guess for fixed and varying parameters, as well as parameter names
BladderParam_03 %BCFunc = BladderFunc_01.m, param, and pn in here

% Starting parameter values
fixed = param; 

% Parameters we are varying (all others fixed)
indexchooseparam = find(chooseparam == 1); 
num_param = length(indexchooseparam); % number of parameters
initialguess = nonzeros(chooseparam.*fixed); % inital guess for parameter we are varying


%% Parameter estimation and generating profile likelihoods 
identifiability_matrix = cell(num_cohorts,3);

for cohort = 1:num_cohorts
    % -------------------------- Specifications ---------------------------
    % data collection times
    t = Timedata;
    t = t(~isnan(t)); 

    % data 
    d = ultrasound_meandata(:,1);
    d = d(~isnan(d));

    % initial (ti) and final (tf) for numerical simulation
    ti = t(1);
    tf = t(end);

    % initial condition
    ic = [squeeze(celltype_meandata(1,cohort,:)); 0];

    % dose for treatment cohort
    ds = cohort_dose(cohort,:);

    % ------------------- Estimate the model parameters -------------------
    [paramests, fval, exitflag] = fminsearch(@(p) ModelCost_ultrasound_04b(BCFunc, ChooseParamFunc_03e(fixed,p,indexchooseparam), t,d,[ti,tf],ic,ds),initialguess,optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000,'TolX',1e-7,'TolFun',1e-7));
    % Note: if exitflag = 1, the function converged to a solution
    %       if exitflag = 0, the number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.MaxFunEvals.
    %       if exitflag = -1, the algorithm was terminated by the output function
    paramests = abs(paramests);
    
    % ------------------- Generate Profile Likelihoods --------------------
    
    % Wrapper function
    costfun = @(p) ModelCost_ultrasound_04b(BCFunc, ChooseParamFunc_03e(fixed,p,indexchooseparam), t,d,[ti,tf],ic,ds);
    threshold = chi2inv(1-alpha,num_param)/2 + fval;
    
    profiles = NaN(numpoints*2,3+num_param,num_param);
    
    for i=1:num_param
        [i]
        % Generate a profile for parameter i, using paramests as the starting
        % value and the fitter to do the parameter estimation:
        profiles(:,:,i) = ProfLike_04c(paramests,i,costfun,profrange,numpoints);
        % each profile has columns: profiled parameter value, resulting
        % cost-function (e.g. RSS) value, any flags from the optimizer, and
        % then columns for each of the other parameter estimates.
    end
    
    identifiability_matrix{cohort,1} = profiles;
    identifiability_matrix{cohort,2} = threshold;
    identifiability_matrix{cohort,3} = paramests;
end


%% Determine 95% confidence intervals for each parameter we're considering
% this helps the plots look nice immediately

confidence_intervals = nan(num_cohorts, length(indexchooseparam), 2);

for cohort = 1:num_cohorts
    profiles  = identifiability_matrix{cohort,1};
    threshold = identifiability_matrix{cohort,2};
    paramests = identifiability_matrix{cohort,3};
    for i = 1:length(indexchooseparam)
        x = profiles(:,1,i);
        y = profiles(:,2,i);

        % Find intersection where profile likelihood crosses threshold
        diffs = y - threshold;
        crossIdx = find(diffs(1:end-1).*diffs(2:end) < 0); 
        if length(crossIdx) < 2 % if less than 2 intersections, parameter is not identifiable, so we don't care about it
            continue
        end
        xIntersect = x(crossIdx) + (threshold - y(crossIdx)).*(x(crossIdx+1) - x(crossIdx)) ./ (y(crossIdx+1) - y(crossIdx));
        if length(xIntersect) > 2 % if more than 2 intersections, grab the closest intersections to paramests
            [~,idx] = sort(abs(xIntersect - paramests(i)));
            xIntersect = xIntersect(idx(1:2)); 
        end
        confidence_intervals(cohort,i,:) = [min(xIntersect), max(xIntersect)];
    end
end


%% Plot Profile Likelihoods of all treatment cohorts

% -------------------------- Figure aesthetics ----------------------------
fs = 20; % Font Size
cohort_names = ["Untreated"; "Gem"; "OT-1"; "Gem+OT-1"];

% treatment cohorts are re-ordered
re_order = [1,3,2,4]; % untx, ot-1, gem, gem/ot-1

% any gem parameters are repositioned to the columns on the right
gem_param = [4,7,11,16];
wgp = ismember(indexchooseparam,gem_param);
reposition_param = [find(wgp==0);find(wgp == 1)]; 

numofrows = ceil(length(indexchooseparam));

% ----------------------------- Plot figure -------------------------------
pfig = figure;
profilefig = tiledlayout(num_cohorts,numofrows,'TileSpacing','compact');
for k = 1:num_cohorts
    cohort = re_order(k); 
    profiles = identifiability_matrix{cohort,1};
    threshold = identifiability_matrix{cohort,2};
    for j = 1:length(indexchooseparam)
        i = reposition_param(j);
        if cohort_dose(cohort,1) ~= 0 || ismember(indexchooseparam(i),gem_param) ~= 1 % when no gem is given, gem parameters are not plotted
            tile_num =((k-1)*length(indexchooseparam)) + j;
            nexttile(tile_num);

            % pull confidence_intervals for plotting
            lb_pl = min(confidence_intervals(:,i,1))*0.9;
            ub_pl = max(confidence_intervals(:,i,2))*1.1;
            if isnan(lb_pl)
                lb_pl = min(profiles(:,1,i));
                ub_pl = max(profiles(:,2,i));
            end

            % plot
            hold on
            pl = plot(profiles(:,1,i),profiles(:,2,i),'k','LineWidth',2); % identifiability curve
            bf = scatter(profiles(numpoints,1,i),profiles(numpoints,2,i),50,'r','filled'); % min point
            ci = plot([lb_pl, ub_pl],[threshold threshold],'r--','LineWidth',2); % red dotted threshhold
            
            % x-axis, y-axis limits, labels and legend
            xlim([lb_pl, ub_pl])
            ylim([0 2*threshold])
            if cohort == num_cohorts
                xlabel(pn{indexchooseparam(i)})
                if i == round(length(indexchooseparam)/2)
                    legend([pl,bf,ci],["Profile Likelihood", "Overall Best Fit", "95% Confidence Interval"],'Location','southoutside')
                end
            else
                set(gca, 'XTickLabel', [])
            end
            if i == 1
                ylabel(cohort_names(cohort));
            else
                set(gca, 'YTickLabel', [])
            end

            grid on
            box on
            set(gca,'LineWidth',1,'FontSize',18,'FontName','Arial')
        end
    end
end

profilefig.YLabel.String = 'Least Squares Error';
profilefig.YLabel.FontSize = fs;
title(profilefig,"Profile Likelihoods",'FontSize',fs);
pfig.Position = [100 200 (length(indexchooseparam)*166)+60 ((num_cohorts+1)*100)];

