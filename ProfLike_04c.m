%% Profile Likelihood Generator 
% Original Source: Marisa Eisenberg (marisae@umich.edu)
% Modified for Bladder Cancer Model in Anderson et al. (2026) --
% "A virtual cohort framework with applications to adoptive cell therapy in bladder cancer"

function profile = ProfLike_04c(params,profindex,costfun,factor,numpoints)
% Definitions
%   params = point in parameter space from which to profile (the parameter estimates)
%   profindex = index of the parameter to be profiled
%   factor = fractional/percent range to profile the parameter over
%   numpoints = half the number of points to sample
%   costfun = this is a cost function of only the parameters (the full vector of them)

% Everythng else (data, ICs, etc.) is fixed for the entire profile so is set outside when 
% costfun is defined. In ProfLike, we'll put a wrapper on costfun that
% will fix the profiled parameter.


% Costfun wrapper
fixcostfun = @(shortparams,profparamval)...
    costfun([shortparams(1:profindex-1);profparamval;shortparams(profindex:end)]); 

% shortparams are all the parameters EXCEPT FOR the profiled one. So, it is
% one parameter less than the usual parameter set. Since we are inputting
% it into the original costfun from identifiability_main_04.m, we need to
% make it back into the full parameter set again, hence the full parameter
% set goes in as: [shortparams(1:profindex-1);profparamval;shortparams(profindex:end)]


% Profile
profrangeDown = linspace(params(profindex), params(profindex)*(1-factor),numpoints)'; % generate (n = numpoints) evenly spaced points between params(profindex) and params(profindex)*(1-factor)
profrangeUp = linspace(params(profindex), params(profindex)*(1+factor),numpoints)';
% split into up and down so we can use last fitted value as starting value for next run
profrange = [profrangeDown profrangeUp];



currfvals = NaN(numpoints*2,1);
currparams = NaN(length(params),numpoints*2);
currflags = NaN(numpoints*2,1);
for i=1:2
    paramstemp = [params(1:profindex-1); params(profindex+1:end)]; % for the rest of the parameters (besides the one being profiled), our initial guess is their previously estimated value
    for j = 1:numpoints % runs through (n = numpoints) different uniformly sampled options for the parameter for the starting points for fminsearch
        [i j] %track progress
        [paramstemp, fvaltemp, flagtemp] = fminsearch(@(p) fixcostfun(p,profrange(j,i)),paramstemp,optimset('MaxFunEvals',5000,'MaxIter',5000)); %re-estimates the other parameters p while keeping the profiled parameter (profrange) fixed 
        currfvals(j+(i-1)*numpoints) = fvaltemp; % objective function values at solutions
        currflags(j+(i-1)*numpoints) = flagtemp; % reason fminsearch terminated
        currparams(:,j+(i-1)*numpoints) = [paramstemp(1:profindex-1);profrange(j,i);paramstemp(profindex:end)]; %storing the profiled value too, so the output parameter values are easy to run the model with
    end
end


currparamsflip = currparams';

profile = [flipud([profrangeDown currfvals(1:numpoints) currflags(1:numpoints) currparamsflip(1:numpoints,:)]);...
    [profrangeUp currfvals(numpoints+1:end) currflags(numpoints+1:end) currparamsflip(numpoints+1:end,:)]];


end