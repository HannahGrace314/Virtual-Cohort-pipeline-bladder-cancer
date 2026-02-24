%% Latin Hypercube Sampling (LHS) for Parameters
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Uses: chooseparam, param, n, and mpr from whatever code is calling this
% script (ParamFitting_03b.m and VirtualCohort_sampling_06.m)

icp = find(chooseparam == 1); % index of chosen parameters 

Psample = (param*ones(1,n))';

rng default
LHSMatrix = lhsdesign(n,length(icp));

for j = 1:length(icp)
    Psample(:,icp(j)) = (mpr(icp(j),2) - mpr(icp(j),1)).*LHSMatrix(:,j) + mpr(icp(j),1);
end

