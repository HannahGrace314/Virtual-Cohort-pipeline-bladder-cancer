%% Numerical simulation of best fit
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Figure code to plot the bladder cancer model's best fit 
% compared to ultrasounds and histology/flow-modified ultrasound data. 
% Uses: BladderParam_03.m and BladderData_shortsummary_02.m

%% Setup

% Cohort aesthetics
cohort_styles = {'-',':', '--', '-.'};
cohort_markers = {'o', 's', '^', 'd'}; 
cohort_names = ["Untreated", "Gem", "OT-1", "Gem+OT-1"];
re_order = [1,3,2,4];

% Cell type aesthetics
 colors = struct( ...
        'cancer', '#22A884', ...
        'tcell', '#4477AA', ...
        'mdsc', '#AA3377', ...
        'gem', [241, 193,0]./255, ... 
        'total', [0.2 0.2 0.2]);

% General aesthetics
tickwidth = 1;
width = 2.5;
fs = 18;
xlim_value = [6,21.5];


%% Figure 5: Numerical simulation of best fit
BladderParam_03
BladderData_shortsummary_02

figure
t = tiledlayout(2, 5, 'TileSpacing','compact', 'Padding','compact');

for cohort = 1:4
    % Load data
    dose = cohort_dose(cohort,:);
    cell_voldata = squeeze(celltype_meandata(:,cohort,:));
    ultrasound_mean = ultrasound_meandata(:,cohort);
    initialcondition = [cell_voldata(1,:), 0]; 

    % Aesthetics
    line_style = cohort_styles{cohort};
    marker_style = cohort_markers{cohort};
    
    % Solve ODE
    [T1,X1] = ode23s(@(t,x) BCFunc(t,x,param), [Timedata(1):0.1:10], initialcondition);
    [T2,X2] = ode23s(@(t,x) BCFunc(t,x,param), [10:0.1:14], [X1(end,1:3), X1(end,4)+dose(1)]);
    [T3,X3] = ode23s(@(t,x) BCFunc(t,x,param), [14:0.1:23], [X2(end,1), X2(end,2)+dose(2), X2(end,3:4)]);
    T = [T1;T2;T3];
    X = [X1;X2;X3];
    
    % Flow: calculate scatter point based on flow ratio of gem-treated to 
    % untreated for gem-treated groups
    if cohort == 2 || cohort == 4
        if dose(1) > 0
            flowdata = [0.413, 0.542];
            [Tgem, Xgem] = ode45(@(t,x) BCFunc(t,x,param), [Timedata(1), 14], initialcondition);
            untxpoint = Xgem(end,2:3);
            day14data = flowdata .* untxpoint;
        end
    end
    
    % --------------------------- Cancer plots ---------------------------- 
    if cohort == 1 || cohort == 3
        ax1 = nexttile(1);
        set(gca,'XTickLabel',{})
        title('Cancer', 'FontSize', fs)
        ylabel('No Gem', 'FontSize', fs-4)
    else
        ax1 = nexttile(1+5);
        ylabel('Gem-treated', 'FontSize', fs-4)
    end

    hold on
    plot(T, X(:,1), 'LineStyle', line_style, 'Color', colors.cancer, 'LineWidth', width)
    scatter(Timedata, cell_voldata(:,1), 'b', 'filled',marker_style)
    
    xlim(xlim_value)
    ylim([0 120])
    set(gca, 'FontSize', fs, 'LineWidth', tickwidth) 
    grid on
    
    % --------------------------- T Cell plot ----------------------------- 
    if cohort == 1 || cohort == 3
        ax2 = nexttile(2);
        set(gca,'XTickLabel',{})
        title('T cells', 'FontSize', fs)
    else
        ax2 = nexttile(2+5);
    end

    hold on
    plot(T1, X1(:,2), 'LineStyle', line_style, 'Color', colors.tcell, 'LineWidth', width)
    plot(T2, X2(:,2), 'LineStyle', line_style, 'Color', colors.tcell, 'LineWidth', width)
    plot(T3, X3(:,2), 'LineStyle', line_style, 'Color', colors.tcell, 'LineWidth', width)
    scatter(Timedata([1,4,end]), cell_voldata([1,4,end],2), 'b', 'filled',marker_style)
    if cohort == 2 || cohort == 4
        scatter(14, day14data(1), 60, 'r', 'filled',marker_style);
    end

    xlim(xlim_value)
    ylim([0 2])
    set(gca, 'FontSize', fs, 'LineWidth', tickwidth)
    grid on
    
    % ---------------------------- MDSC plot ------------------------------ 
    if cohort == 1 || cohort == 3
        ax3 = nexttile(3);
        set(gca,'XTickLabel',{})
        title('MDSCs', 'FontSize', fs)
    else
        ax3 = nexttile(3+5);
    end

    hold on
    plot(T, X(:,3), 'LineStyle', line_style, 'Color', colors.mdsc, 'LineWidth', width)
    scatter(Timedata([1,4,end]), cell_voldata([1,4,end],3), 'b', 'filled',marker_style)
    if cohort == 2 || cohort == 4
        sc_flow = scatter(14, day14data(2), 60, 'r', 'filled',marker_style);
        sc_data = scatter(Timedata([4,end]), cell_voldata([4,end],3), 40, 'b', 'filled',marker_style);
    else
        sc_data = scatter(Timedata([4,end]), cell_voldata([4,end],3), 40, 'b', 'filled',marker_style);
    end

    xlim(xlim_value)
    ylim([0 10])
    set(gca, 'FontSize', fs, 'LineWidth', tickwidth)
    grid on
    
    % ------------------------- Total Tumor plot -------------------------- 
    if cohort == 1 || cohort == 3
        ax4 = nexttile(4);
        set(gca,'XTickLabel',{})
        title('Total Tumor', 'FontSize', fs)
    else
        ax4 = nexttile(4+5);
    end

    hold on
    plot(T1, sum(X1(:,1:3),2), 'LineStyle', line_style, 'Color', colors.total, 'LineWidth', width);
    plot(T2, sum(X2(:,1:3),2), 'LineStyle', line_style, 'Color', colors.total, 'LineWidth', width)
    plot(T3, sum(X3(:,1:3),2), 'LineStyle', line_style, 'Color', colors.total, 'LineWidth', width)
    % Ultrasound data 
    h_ultrasound(cohort) = scatter(Timedata, ultrasound_mean, 'k', 'filled',marker_style);
    % For Cohort Legend
    allcohorts_leg(cohort) = plot(nan, nan, sprintf('k%s%s',line_style,marker_style), 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k','LineWidth', width/2);
    
    xlim(xlim_value)
    ylim([0 120])
    set(gca, 'FontSize', fs, 'LineWidth', tickwidth)
    grid on
    
    % ----------------------------- Gem plot ------------------------------ 
    if cohort == 1 || cohort == 3
        ax5 = nexttile(5);
        set(gca,'XTickLabel',{})
        title('Gem', 'FontSize', fs)
    else
        ax5 = nexttile(5+5);
    end

    hold on
    plot(T1, X1(:,4), 'LineStyle', line_style, 'Color', colors.gem, 'LineWidth', width)
    plot(T2, X2(:,4), 'LineStyle', line_style, 'Color', colors.gem, 'LineWidth', width)
    plot(T3, X3(:,4), 'LineStyle', line_style, 'Color', colors.gem, 'LineWidth', width)

    set(gca, 'FontSize', fs, 'LineWidth', tickwidth)
    grid on
    ylim([0 4e4])
    set(gca,'YTick', [1,2,3,4]*10^4,'YTickLabel',{'1','2','3','4e4'})
    ylabel(ax5,'Concentration (\muM)','FontSize',fs)
    xlim(xlim_value)
end

% Treatment groups legend
leg7 = legend(ax4, allcohorts_leg(re_order), cohort_names(re_order), 'FontSize', fs-2);
leg7.Layout.Tile = 'east';
leg7.Title.String = 'Treatment Groups';

% Data legend
leg5 = legend(ax5, [sc_data, sc_flow, h_ultrasound(4)], {'Histology-modified Ultrasound', 'Flow Cytometry','Ultrasound'}, 'FontSize', fs-2);
leg5.Layout.Tile = 'east';
leg5.Title.String = 'Data Types';

% Final titles and labels
set(gcf, 'Position', [100, 200, 1200, 400]);
xlabel(t,'Time (days)','FontSize',fs+2)
ylabel(t,'Volume (mm^3)','FontSize',fs+2)
