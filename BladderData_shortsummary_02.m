%% Short Summary of Bladder Cancer Data 
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Bladder cancer data for the 4 experimental treatment cohorts 
% (untreated, Gem, OT-1, Gem+OT-1) includes:
% (1) data collection types (Timedata),
% (2) volume data for each cell type (celltype_meandata and celltype_stddata),
% (3) longitudinal ultrasound data for total tumor size (ultrasound_meandata and ultrasound_stddata),
% (4) treatment injections (dose_cohort), and
% (5) initial percentage of cells at day 6 from linear interpolation 
% (interp_percent(1,:,1)) in BladderData_analysis_02c.m.

%% Data for an average mouse in each treatment cohort
Timedata = [6; 9; 13; 16; 20];
celltype_meandata = zeros(length(Timedata),4,3);
celltype_stddata = zeros(length(Timedata),4,3);

exceltitle = ["cancer percent (C = tumor-T-M)","treatment T cell percent","treatment MDSC percent"];

% Ultrasound data (mm^3) modified by Histology and Flow for each cell type
for celltype = 1:3
    celltype_meandata(:,:,celltype) = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_ODE_VC_GitHub/histology_flow_modified_ultrasound_data.xlsx','Sheet',exceltitle(celltype),'Range','B48:E52'));
    celltype_stddata(:,:,celltype) = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_ODE_VC_GitHub/histology_flow_modified_ultrasound_data.xlsx','Sheet',exceltitle(celltype),'Range','B59:E63'));
    % untreated, gem mono, ot1 mono, gem+ot1 combo--each have a column for the 5 time points
end

% Ultrasound data (mm^3)
ultrasound_meandata = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_ODE_VC_GitHub/histology_flow_modified_ultrasound_data.xlsx','Sheet',exceltitle(1),'Range','B26:E30'));
ultrasound_stddata = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_ODE_VC_GitHub/histology_flow_modified_ultrasound_data.xlsx','Sheet',exceltitle(1),'Range','B37:E41'));
% untreated, gem mono, ot1 mono, gem+ot1 combo--each have a column for the 5 time points


%% Treatment injection based on treatment cohort

prot1 = 0.75; % percent remaining after OT-1 injection (based on human data - private communication)

cohort_dose = [0,     0;            % cohort 1 - untreated
               3.8e4, 0;            % cohort 2 - gem mono (gem dose converted into uM)
               0,     1.426*prot1;  % cohort 3 - ot1 mono (1.426 mm^3 = 5e6 OT-1 cells converted to mm^3 based on T cell volume estimate from day 17 histology)
               3.8e4, 1.426*prot1]; % cohort 4 - gem+ot1 combo


%% Initial percentage of cells at day 6 from linear interpolation

initialpercent = [0.9710,0.0026,0.264];
