%% Bladder Cancer Data Analysis
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Script for Anderson et al. (2026) which 
% (1) re-groups ultrasound data and plots it (Figure(1))
% (2) finds the maximum number of standard deviations needed to test for virtual cohort sampling,
% (3) plots histology data at day 17 and 23 (Figure(3)),
% (4) performs a linear interpolation of percentage data from histology 
% and flow cytometry and plots it (Figure(4)), and 
% (5) captures cell percentages in terms of the ultrasound data collection 
% times for enhanced model fitting (as per structural identifiability result).

%% Figure aesthetics

fs = 15; % FontSize

% Cohort:
cohort_names = ["Untreated","Gem", "OT-1","Gem+OT-1"];
cohort_markers = {'o', 's', '^', 'd'}; 
cohort_styles = {'-',':', '--', '-.'};
re_order = [1,3,2,4]; % for re-ordering plots of cohorts

% Cell:
cell_colors = ["#22A884", "#4477AA", "#AA3377"]; % Cancer, T cells, MDSCs
cell_names = ["Cancer"; "T cell"; "MDSC"];


%% Volume data in terms of mm^3 

BladderData_allmice_02b


%% Grouping data according to different treatment groups
% we only start grouping data separately once treatment starts
num_times = length(Timedata);

groupData = cell(num_times, 4);
meandata = zeros(num_times,4);
stddata = zeros(num_times,4);
mindata = zeros(num_times,4);
maxdata = zeros(num_times,4);

for day = 1:num_times
    for cohort = 1:4
        if day <= 2  % Days 6 and 9: all groups share same data since all untreated at this point
            vec = squeeze(exp_individual(:,day,:));
        elseif day == 3  % Day 13: Gem treatment at day 10, so this data split between Gem-treated and Gem-untreated
            idx = [2,4] - mod(cohort,2);
            vec = [exp_individual(idx,day,:), exp_individual(idx,day,:)];
        else  % Days 16, 20: OT-1 treatment at day 14, so all four cohorts (untreated, Gem, OT-1, Gem+OT-1) are distinct now
            vec = exp_individual(cohort,day,:);
        end

        groupData{day,cohort} = vec(~isnan(vec));

        % Calculate the mean and standard deviation for all re-grouped cohorts
        meandata(day,cohort) = mean(groupData{day,cohort});
        stddata(day,cohort) = std(groupData{day,cohort});
        mindata(day,cohort) = min(groupData{day,cohort}); 
        maxdata(day,cohort) = max(groupData{day,cohort}); 
    end
end


%% Figure (unused): Tumor volume and errorbars for the four treatment cohorts
% we only start grouping data separately once treatment starts
plot_color = [0.2,0.2,0.2];

% Plot figure
figure
hold on
tumorline = gobjects(1, 4);
for cohort = 4:-1:1 
    style = cohort_styles{cohort}; 
    color = plot_color*(-cohort+4);

    % Plot lines
    tumorline(cohort) = plot(Timedata, meandata(1:end, cohort), ...
        'LineWidth', 5, 'Color', color, 'LineStyle', style);  

    % Plot errorbars

    % Untreated
    if cohort == 1
        errorbar(Timedata, meandata(1:end, cohort), stddata(1:end, cohort), cohort_markers{cohort}, ...
            'LineWidth', 2, 'Color', color);  
    end

    % Gem treated
    if cohort == 2
        errorbar(Timedata(3:end), meandata(3:end, cohort), stddata(3:end, cohort), cohort_markers{cohort}, ...
            'LineWidth', 2, 'Color', color);  
    end

    % OT-1 / Gem+OT-1 treated
    if cohort == 3 || cohort == 4
        errorbar(Timedata(4:end), meandata(4:end, cohort), stddata(4:end, cohort), cohort_markers{cohort}, ...
            'LineWidth', 2, 'Color', color); 
    end
    % Note: errorbars are calculated so that until the treatment regimens
    % start to differ, mice are grouped together (regardless of if they
    % will eventually be in a different cohort)
    % Example: at day 6, no mice have been treated, so the data from all
    % 55 mice creates the errorbar
end

legend(tumorline(re_order), cohort_names(re_order), 'Location', 'Northwest');

ylim([0, 230]);
ylabel('Tumor Volume (mm^3)');
xlabel('Days');
set(gca,'XTick', Timedata)
grid on;

sgtitle("Ultrasound: Total Tumor Volume", 'FontSize', fs,'FontWeight','bold');
set(gca, 'FontSize', fs);


%% Determining the maximum number of stds for virtual cohort sampling
A = maxdata - meandata;
B = meandata - mindata;
max_stds = max(ceil(max(A,B)./stddata),[],'all');


%% Load and organize histology data based on treatment cohort

% Load data 
excelpercentdata = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_ODE_VC_GitHub/histology_flow_data.xlsx','Sheet','collective data histology','Range','B3:J6'));

percentdata_day17 = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_ODE_VC_GitHub/histology_flow_data.xlsx','Sheet','day 17','Range','F3:H14'));

percentdata_day23 = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_ODE_VC_GitHub/histology_flow_data.xlsx','Sheet','day 23','Range','F3:H14'));

% Organize data (percent data in mean_histologydata and all_histologydata out of 1 (instead of 100))
percenttimedata = [0;17;23];
mean_histologydata = nan(3,3,4); 
all_histologydata = cell(2,3,4); 
for cohort = 1:4
    mean_histologydata(:,:,cohort) = reshape(excelpercentdata(cohort,:),[3,length(excelpercentdata(cohort,:))/3])';
    % mean_histologydata(:,:,cohort) is a 3x3 matrix: 
    % column 1 includes percentages of cancer cells in the tumor at day 0, 17, and 23
    % column 2 includes percentages of T cells "
    % column 3 includes percentages of MDSCs "

    for celltype = 1:3
        day17 = percentdata_day17(((cohort*3)-2):(cohort*3),celltype);
        day23 = percentdata_day23(((cohort*3)-2):(cohort*3),celltype);
        all_histologydata{1,celltype,cohort} = day17(~isnan(day17))';
        all_histologydata{2,celltype,cohort} = day23(~isnan(day23))';
    end
    % all_histologydata is a 2x3x4 cell representing 
    % day (17 or 23) x cell type (cancer, T cell, MDSC) x treatment cohort 
    % (untreated, Gem, OT-1, Gem+OT-1)
end


%% Figure 3: Percentage of MDSCs and T cells in Tumor from Histology 

% Grab histology max and min at each day (17 or 23) and treatment cohort
[errorbarmin, errorbarmax] = cellfun(@(x) bounds(x * 100), all_histologydata);

figure(3)
tiledlayout(2,4,"TileSpacing","compact","Padding","tight");
for j = 1:4
    cohort = re_order(j);

    % CD8+ T cells
    nexttile(j);
    grid on
    hold on
    h(2) = errorbar([1,3], mean_histologydata(2:3,2,cohort)*100,...
        mean_histologydata(2:3,2,cohort)*100 - errorbarmin(:,2,cohort),...
        errorbarmax(:,2,cohort) - mean_histologydata(2:3,2,cohort)*100, 'o',"LineWidth",2, "Color",cell_colors(3));
    h(1) = errorbar([1,3], mean_histologydata(2:3,2,cohort)*100,...
        mean_histologydata(2:3,2,cohort)*100 - errorbarmin(:,2,cohort),...
        errorbarmax(:,2,cohort) - mean_histologydata(2:3,2,cohort)*100, 'o',"LineWidth",2, "Color",cell_colors(2));
    xticks([1 3]);
    yticks([0 0.5 1 1.5 2])
    xticklabels({'Day 17','Day 23'});
    xlim([0 4])
    ylim([0 2])
    if cohort == 1
        ylabel('% CD8 in tumor');
        legend(h, ["CD8+ T cells","Ly6G+ MDSCs"],'Location','northeast')
    end
    title(cohort_names(cohort),'FontWeight','normal');

    % Ly6G+ MDSCs
    nexttile(j+4);
    grid on
    hold on
    errorbar([1,3], mean_histologydata(2:3,3,cohort)*100,...
        mean_histologydata(2:3,3,cohort)*100 - errorbarmin(:,3,cohort),...
        errorbarmax(:,3,cohort) - mean_histologydata(2:3,3,cohort)*100, 'o',"LineWidth",2, "Color",cell_colors(3));
    xticks([1 3]);
    xticklabels({'Day 17','Day 23'});
    xlim([0 4])
    ylim([0 15])
    if cohort == 1
        ylabel('% Ly6G in tumor');
    end
end

sgtitle('Percentage of MDSCs and T cells in Tumor from Histology','FontWeight','bold')


%% Straight line interpolation based on T cell and MDSC histology and flow data 
flowdata = nan(1,3); % cancer, T cells, MDSCs - ratio of (gem treated / untreated) at day 14
flowdata(2:3) = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_ODE_VC_GitHub/histology_flow_data.xlsx','Sheet','day 14','Range','B7:C7'));
flow_data_time = 14;
treat_days = [10;14];
wantedtimes = sort(unique([Timedata;percenttimedata;treat_days;flow_data_time]));

linearinterp = zeros(length(wantedtimes),3,4);

for i = 2:3 % T cells and MDSCs
    % Untreated
    linearinterp(:,i,1) = interp1(percenttimedata',mean_histologydata(:,i,1),wantedtimes');

    % Gem mono - use day 6, 9, 10 from untreated
    linearinterp(:,i,2)= interp1([wantedtimes(1:4);flow_data_time;percenttimedata(2:3)]',...
        [linearinterp(1:4,i,1);linearinterp(6,i,1)*flowdata(i);mean_histologydata(2:3,i,2)],wantedtimes');

    % OT-1 mono - use day 6, 9, 10, 13, 14 from untreated
    linearinterp(:,i,3) = interp1([wantedtimes(1:6);percenttimedata(2:3)]',...
        [linearinterp(1:6,i,1);mean_histologydata(2:3,i,3)],wantedtimes');

   % Gem+OT-1 combo - use day 6, 9, 10, 13, 14 from Gem mono
    linearinterp(:,i,4) = interp1([wantedtimes(1:6);percenttimedata(2:3)]',...
        [linearinterp(1:6,i,2);mean_histologydata(2:3,i,4)],wantedtimes');
end

% Cancer: C = tumor - T - M
for cohort = 1:4
    linearinterp(:,1,cohort) = 1-sum(linearinterp(:,2:3,cohort),2);
end


%% Figure 4: Linear Interpolation of Cell Percentages from Histology and Flow 

figure(4)
linearinterp_fig = tiledlayout(1,3, 'TileSpacing', 'Compact', 'Padding', 'Compact'); 
for i = 1:3 % cell types
    nexttile(i)
    hold on
    treat_group = gobjects(1,4);
    color = cell_colors(i);
   
    for cohort = 1:4 
        style = cohort_styles{cohort};

        plot(wantedtimes, linearinterp(:,i,cohort) * 100,'LineWidth', 3, 'Color', color, 'LineStyle', style)

        % assumption at day 0
         assump_leg = scatter(wantedtimes(1),linearinterp(1,i,cohort)*100, ...
            60, 'MarkerEdgeColor',[34, 139, 34]./255,'MarkerFaceColor',[34, 139, 34]./255, 'MarkerFaceAlpha', 0.8);

         if i == 1 % Cancer
            % for Legend 1: cohorts
            treat_group(cohort) = plot(nan, nan, ...
               'LineWidth', 3, ...
               'Color', 'k', ...
               'LineStyle', style); 

        else % T cell and MDSCs
            % T cells & MDSCs histology scatter points
            hist_leg = scatter(percenttimedata(2:3), mean_histologydata(2:3,i,cohort)*100, ...
                60, 'o','b', 'filled', 'MarkerFaceAlpha', 0.8); 

            % T cells & MDSC scatter point affected by flow cytometry
            flow_leg = scatter(flow_data_time, linearinterp(wantedtimes == flow_data_time,i,2)*100, ...
                60, 'o','r', 'filled', 'MarkerFaceAlpha', 0.8); 
         end
          
    end

    % Legend 1: cohorts
    if i == 1
        legend(treat_group(re_order),cohort_names(re_order),'Location', 'southwest', 'FontSize', fs-2)
    end
    
    % Legend 2: data type
    if i == 3
        legend([hist_leg,flow_leg,assump_leg],["Histology","Flow","Assumption"],'Location', 'northwest', 'FontSize', fs-2)
    end

    set(gca, 'FontSize', fs, 'LineWidth', 1.5);
    ylabel(join(["% " cell_names(i)]), 'FontSize', fs);
    grid on
end

set(gcf,'Position',[0 300 900 300],'PaperPositionMode','auto');
xlabel(linearinterp_fig,'Days', 'FontSize', fs);
sgtitle('Linear Interpolation of Cell Percentages', 'FontSize', fs, 'FontWeight', 'bold');


%% Percentages from interpolation for ultrasound data collection times
interp_percent = linearinterp(ismember(wantedtimes,Timedata),:,:); % ultrasound times x celltype x cohort