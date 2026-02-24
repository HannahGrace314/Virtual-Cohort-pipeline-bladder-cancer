%% Numerical simulation figure of virtual cohort
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Using choice_std, choose how many standard deviations from 
% mean within which to keep virtual mice and then code pulls those 
% virtual mice, simulates these mice for any treatment cohorts that were 
% left out of sampling, and then plots the numerical simulations of the
% virtual cohort alongside data for each treatment cohort.

% Uses: BladderData_shortsummary_02.m, BladderFunc_01.m,
% and Psample, accepted_withinstd_all, collect_numsim, initialcondition_n
% and timerange from output of VirtualCohort_sampling_06.m

%% Parameter sets accepted based on number of standard deviations from the mean

choice_std = 2;
% can do up to 4

vc_param = Psample(accepted_withinstd_all(:,choice_std),:);
idx_accepted = find(accepted_withinstd_all(:,choice_std)==1);

num_cohorts = 4;


%% ODE evaluation 
BladderData_shortsummary_02 % has cohort_dose, celltype_meandata (for initial condition), ultrasound_meandata, and Timedata

BCFunc = @BladderFunc_01;
n = length(vc_param);
vc_initialcondition = initialcondition_n(accepted_withinstd_all(:,choice_std),:);
vc_numsim = nan(num_cohorts,4,length(timerange), n); % cohort x celltype x time x mouse

for cohort = 1:num_cohorts
    dose = cohort_dose(cohort,:);
    if cohort <= num_samplecohort % pulling simulations of accepted mice from virtual cohort sampling
        vc_numsim(cohort,:,:,:) = collect_numsim(cohort,:,:,idx_accepted);
    else  % evaluating experimental cohorts not involved in virtual cohort sampling
        parfor mouse = 1:n
            pull_numsim = Evaluate_NumSim_03g(BCFunc,vc_param(mouse,:),[Timedata(1),23],vc_initialcondition(mouse,:),dose);
            vc_numsim(cohort,:,:,mouse) = [pull_numsim{1}, sum(pull_numsim{1},2)]';
        end
    end
end


%% Capturing flow data point (compared to untreated simulation using param set from parameter estimation)

flowdata = [0.413, 0.543];
initialcondition = [squeeze(celltype_meandata(1,1,:)); 0];
[Tgem, Xgem] = ode45(@(t,x) BCFunc(t,x,param), [Timedata(1), 14], initialcondition);
untxpoint = Xgem(end,2:3);
day14data = flowdata .* untxpoint;


%% Figure: Comparison of all cell types to all treatment groups

fs = 15;

% Cell type aesthetics
colors = [34/255, 168/255, 132/255; % cancer
         68/255, 119/255, 170/255; % T cell
         170/255, 51/255, 119/255; % MDSC
         0.5 0.5 0.5]; % total tumor
yaxis_labels = ["Cancer", "T cells", "MDSCs", "Total Tumor"];
ylim_range = [0 300; 0 5; 0 25; 0 300];

% Cohort aesthetics
cohort_names = ["Untreated","Gem","OT-1","Gem+OT-1"];
re_order = [1,3,2,4];

figure
allgroups_virtualcohort = tiledlayout(4,num_cohorts,'TileSpacing','compact'); 
for cohort = 1:num_cohorts
    j = re_order(cohort);

    % Load data
    cell_voldata = squeeze(celltype_meandata(:,cohort,:));
    ultrasound_mean = ultrasound_meandata(:,cohort);
    ultrasound_std = ultrasound_stddata(:,cohort);

    for celltype = 1:4
        nexttile(j+(4*(celltype-1)))
        hold on

        % Simulations
        for mouse = 1:n
            plot(timerange, squeeze(vc_numsim(cohort,celltype,:,mouse)),'Color',colors(celltype,:)*(mouse/n));
        end
        
        % Data
        if celltype == 1
            scatter(Timedata, cell_voldata(:,1), 'b', 'filled')
            title(cohort_names(cohort), 'FontWeight', 'bold','FontSize',fs-2)
        end
        if celltype == 2 || celltype == 3
            scatter(Timedata([4,end]), cell_voldata([4,end],celltype), 'b', 'filled');
            if cohort == 2 || cohort == 4
                scatter(14,day14data(celltype-1),'r','filled'); % for gem/untx only cohort
            end
        end
        if celltype == 4
            scatter(Timedata(1), ultrasound_mean(1), 'k','filled');
            errorbar(Timedata, ultrasound_mean, choice_std*ultrasound_std, 'k', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'k','LineWidth',2);
            set(gca,'XTick', Timedata,'XTickLabel',Timedata,'FontSize',fs-4)
        else
            set(gca,'XTick', Timedata,'XTickLabel',{},'FontSize',fs-4)
        end

        if cohort == 1
            ylabel(yaxis_labels(celltype))
        else
            set(gca,'YTickLabel',{})
        end

        ylim(ylim_range(celltype,:))
        xlim([6,21.5])
        grid on
    end
end

title(allgroups_virtualcohort,'Virtual Murine Cohort','FontSize',fs+2)
xlabel(allgroups_virtualcohort,'Time (days)','FontSize',fs)
ylabel(allgroups_virtualcohort,'Volume (mm^3)','FontSize',fs)
set(gcf,'Position',[0 300 1000 800],'PaperPositionMode','auto');




