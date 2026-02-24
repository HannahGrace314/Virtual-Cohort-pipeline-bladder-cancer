%% Dynamics Validation
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Finds the virtual mouse that best replicates the dynamics 
% of each original experimental mouse.
% Uses: Bladder_Data_allmice_02b.m, BladderFunc_01.m, 
% ModelCost_ultrasound_04b.m, and the output from virtualcohort_numsim_06b.m
% (namely, variables vc_param, vc_initialcondition, and num_cohorts)

%% Pull ODE and ultrasound data

% ODE
BCFunc = @BladderFunc_01;

% Data
BladderData_allmice_02b % ultrasound data for all mice
BladderData_shortsummary_02 % has cohort doses


%% Find closest virtual mouse to each experimental mouse in terms of least squares error
mouse_errors = nan(length(vc_param),num_cohorts,15);
vmouse_numsim = cell(num_cohorts,15);

% Choose a experimental mouse
for cohort = 1:num_cohorts
    for emouse = 1:15
        data = rmmissing(exp_individual(cohort,:,emouse))';
        times = rmmissing(times_all(cohort,:,emouse));
    
        if ~isempty(times)
            % ------ Compute error of each virtual mice to experimental -------
            for vmouse = 1:length(vc_param)
                param = vc_param(vmouse,:);
                mouse_errors(vmouse,cohort,emouse) = ModelCost_ultrasound_04b(BCFunc,param, times,data,[6,23],vc_initialcondition(vmouse,:),cohort_dose(cohort,:)); 
            end
    
            % --- Simulate closest virtual mouse to each experimental mouse ---
            vmouse = find(mouse_errors(:,cohort,emouse) == min(mouse_errors(:,cohort,emouse)));
            param = vc_param(vmouse,:);
            numsim = Evaluate_NumSim_03g(BCFunc,param,[6,23],vc_initialcondition(vmouse,:),cohort_dose(cohort,:));
            vmouse_numsim{cohort,emouse} = {sum(numsim{1},2),numsim{2}}; % {total tumor, time}
        end
    end
end


%% Plot closest virtual mouse to each experimental mouse 
   
% --------------------------- Figure aesthetics ---------------------------
tickwidth = 0.8;
width = 2;
fs = 17;
cohort_names = ["Untreated","Gem", "OT-1", "Gem+OT-1"];
plotcolor = [34/255,  168/255, 132/255; % green
             68/255,  119/255, 170/255; % blue
             170/255,  51/255, 119/255; % pink
             241/255, 193/255,       0; % yellow
             0.5 0.5 0.5;               % dark grey
             147/255, 145/255,       1; % purple
             238/255, 119/255,  51/255; % orange
             0 0.95 0.95;               % cyan
             0.2, 0.9, 0.7;             % bright green
             0.4, 0.7, 1.0;             % bright blue
             1.0, 0.2, 0.6;             % neon pink
             1.0, 0.8, 0.0;             % neon yellow
             0.8, 0.8, 0.8;             % light grey
             1.0, 0.6, 0.2;             % neon orange
             0.7, 0.7, 1.0;             % bright purple
             0 0.6 0.6];                % teal

% -------------------- Plot dynamics comparison figure --------------------
figure
numsimplot = tiledlayout(2,num_cohorts,'TileSpacing','compact');

for cohort = 1:num_cohorts
    % ---------------------------- Virtual mice ---------------------------- 
    nexttile(cohort);
    hold on
    for emouse = 1:15
        if ~isempty(vmouse_numsim{cohort,emouse})
            x = vmouse_numsim{cohort,emouse}{1}(ismember(vmouse_numsim{cohort,emouse}{2}, times_all(cohort,:,emouse)), :); 
            t = nonzeros(ismember(vmouse_numsim{cohort,emouse}{2}, times_all(cohort,:,emouse)).*vmouse_numsim{cohort,emouse}{2}); 
            plot(t,x,'Color',plotcolor(emouse,:),'LineWidth',width)
        end
    end
    ylim([0 250])
    xlim([9 20])
    set(gca,'FontSize',fs,'LineWidth',tickwidth)
    title(cohort_names(cohort))
    xticklabels([])
    if cohort>1
        yticklabels([])
    else
        ylabel('Virtual')
    end
   
    % ------------------------- Experimental mice ------------------------- 
    nexttile(cohort+num_cohorts);
    hold on
    for emouse = 1:15
        if ~isempty(times_all(cohort,:,emouse))
            plot(times_all(cohort,:,emouse), exp_individual(cohort,:,emouse),'Color',plotcolor(emouse,:),'LineWidth',width); % NEED TO FIX times_all here
        end
    end
    ylim([0 250])
    xlim([9 20])
    set(gca,'FontSize',fs,'LineWidth',tickwidth)
    if cohort>1
        yticklabels([])
    else
        ylabel('Experimental')
    end
end

xlabel(numsimplot, 't (days)','FontSize',fs+2)
ylabel(numsimplot, 'Total tumor (mm^3)','FontSize',fs+2)
set(gcf,'Position',[0 300 1000 540],'PaperPositionMode','auto');
title(numsimplot,'Closest virtual mouse to each experimental mouse','FontWeight', 'bold','FontSize',fs+2)

