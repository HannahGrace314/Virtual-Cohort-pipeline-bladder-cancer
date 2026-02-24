%% Virtual Cohort Distribution Validation
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Script determines the summary statistics of the virtual and
% experimental cohorts, plots the violin plot comparing these, and runs the
% Kolmogorov-Smirnov and Student's t test to determine if there is a
% significant difference between these.

% Uses: BladderData_allmice_02b.m and output (vc_numsim) from virtualcohort_numsim_06b.m

num_cohorts = 4;

%% Summary Statistics of Virtual Cohort

[~, idx] = ismember(Timedata, timerange);

% Virtual Cohort Stats
vc.stats = summary_stats(squeeze(vc_numsim(:,4,idx,:))); % input: cohorts x times x number of mice
% index 4 to grab total tumor size


%% Summary Statistics of Experimental Cohort

BladderData_allmice_02b

% Experimental Cohort Stats
data.stats = summary_stats(exp_individual); % input: cohorts x times x number of mice


%% Violin plots comparing virtual and experimental cohorts
cohort_names = ["Untreated","Gem","OT-1","Gem+OT-1"];
fs = 19;
vc_color = [241, 193,0]./255; 
ec_color = [0,0,0];

figure
vp = tiledlayout(2,ceil(num_cohorts/2),"TileSpacing","compact");
for cohort = 1:num_cohorts
    nexttile(cohort)
    violinplot(squeeze(exp_individual(cohort,:,:))',num2cell(Timedata),'ViolinColor',ec_color,'ShowWhiskers',false)
    violinplot(squeeze(vc_numsim(cohort,4,idx,:))',num2cell(Timedata),'ViolinColor',vc_color,'ShowData',false,'ShowWhiskers',false)
    if cohort == 1
        vc = fill(nan(4,1),nan(4,1), vc_color,'FaceAlpha', 0.3, 'EdgeColor', vc_color);
        ec = fill(nan(4,1),nan(4,1), ec_color,'FaceAlpha', 0.3, 'EdgeColor', ec_color);
        data = scatter(nan,nan,[],ec_color,'filled');
        legend([vc,ec,data],["Virtual Cohort","Experimental Cohort", "Data Point"],'Location','northwest')
    end
    subtitle(cohort_names(cohort),'FontSize',fs,'FontWeight','bold');
    set(gca,'FontSize',fs-3,'ylim',[0, 300],'LineWidth',0.75)
    box off
end
xlabel(vp,'Days','FontSize',fs)
ylabel(vp,'Volume (mm^3)','FontSize',fs)


%% Kolmogorov-Smirnov (compares distributions) and Student's t test (compares means)
% for statistically comparing virtual and experimental cohorts
h_ks = nan(num_cohorts, length(Timedata));
p_ks = nan(num_cohorts, length(Timedata));
h_ttest = nan(num_cohorts, length(Timedata));
p_ttest = nan(num_cohorts, length(Timedata));

for cohort = 1:num_cohorts
    virtual_data = squeeze(vc_numsim(cohort,4,idx,:));
    exp_data = squeeze(exp_individual(cohort,:,:));
    exp_data(all(isnan(exp_data),2),:) = []; % removes rows with no mice data

    for day = 1:length(Timedata)
        data1 = virtual_data(day,:); 
        data2 = exp_data(day,~isnan(exp_data(day,:)));

        % Compare full distributions (Kolmogorov-Smirnov)
        [h_ks(cohort,day),p_ks(cohort,day)] = kstest2(data1, data2); 
        % Want p_ks >= 0.05 (so no significant difference found), and 
        % h_ks = 0 (meaning fail to reject the null hypothesis, ie the data comes from similar population)

        % Compare means (Student's t-test)
        [h_ttest(cohort,day),p_ttest(cohort,day)] = ttest2(data1, data2);
        % again, p_ttest >= 0.05 and h_ttest = 0
    end

end


%% Calculate Stats
function stats = summary_stats(input)
% input is a matrix of num_cohorts x tsize (num of timepts) x number of mice
    num_cohorts = length(input(:,1,1));
    tsize = length(input(1,:,1));
    
    stats.mean   = zeros(num_cohorts, tsize);
    stats.median = zeros(num_cohorts, tsize);
    stats.max    = zeros(num_cohorts, tsize);
    stats.min    = zeros(num_cohorts, tsize);
    stats.iqr    = zeros(num_cohorts, tsize);
    stats.std    = zeros(num_cohorts, tsize);

    for c = 1:num_cohorts
        for t = 1:tsize
            grab_data = squeeze(input(c,t,:));
            grab_data = grab_data(~isnan(grab_data)); % remove NaNs if any
            
            stats.mean(c,t)   = mean(grab_data);
            stats.median(c,t) = median(grab_data);
            stats.max(c,t)    = max(grab_data);
            stats.min(c,t)    = min(grab_data);
            stats.iqr(c,t)    = iqr(grab_data);
            stats.std(c,t)    = std(grab_data);
        end
    end
end
