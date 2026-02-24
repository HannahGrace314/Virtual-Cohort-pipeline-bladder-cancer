%% Relative Error Cost Function for model compared to histology/flow data
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Function to calculate relative error of the bladder cancer 
% model compared to histology and flow cytometry data.

function rel_error = ModelCost_histologyflow_03d(BCFunc, params, timerange, data, initialcondition, dose)
params = abs(params);

if params(4)<= min(params(7),params(11)) %kGC <= min{kGT, kGM}
    times = squeeze(data(:,1,:));
    unique_times = sort(unique(times(~isnan(times))));
    voldata = squeeze(data(:,2,:));
    num_datavariables = length(voldata(1,:));

    numsim = Evaluate_NumSim_03g(BCFunc,params,timerange,initialcondition,dose);
    X = numsim{1};
    T = numsim{2};

    % Update num sim values 
    x = X(ismember(T, unique_times), :);
    t = nonzeros(ismember(T, unique_times).*T);

    % Flow data at day 14
    if dose(1)>0 % gem treated 
        flowdata = [NaN, 0.413,0.542]; % ratios that gem decreases CD8 T cells and MDSCs at day 14 compared to non-gem treated (Bazargan et al. 2023, Fig 4BC)
        [~,X_untx] = ode45(@(t,x) BCFunc(t,x,params), [timerange(1),14], initialcondition);
        day14data = flowdata.*X_untx(end,1:num_datavariables);
    end

    % Relative error
    rel_error = 0;
    for v = 1:num_datavariables
        if dose(1) > 0 && v~=1 % when treating with gemcitabine, includes flow data for T and M
            rel_error = rel_error +sum((abs(voldata(t == times(:,v),v) - x(t == times(:,v),v)))./voldata(t == times(:,v),v)) + abs(day14data(v)-X(T==14,v))./day14data(v);
        else
            rel_error = rel_error +sum((abs(voldata(t == times(:,v),v) - x(t == times(:,v),v)))./voldata(t == times(:,v),v));
        end
    end
   
else 
    rel_error = 100;
end
