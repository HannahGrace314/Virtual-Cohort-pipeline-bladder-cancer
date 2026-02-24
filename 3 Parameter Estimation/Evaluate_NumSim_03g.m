%% ODE evaluation function 
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Simulates the bladder cancer model (BCFunc) with a
% gemcitabine dose at day 10 and an OT-1 dose at day 14.

function numsim = Evaluate_NumSim_03g(BCFunc,params,timerange,initialcondition,dose)
    params = abs(params);

    % EVALUATE ODE:
    [T1,X1] = ode45(@(t,x) BCFunc(t,x,params), timerange(1):1:10, initialcondition);

    [T2,X2] = ode45(@(t,x) BCFunc(t,x,params), 10:1:14, [X1(end,1:3),X1(end,4)+dose(1)]);

    [T3,X3] = ode45(@(t,x) BCFunc(t,x,params), 14:1:timerange(end), [X2(end,1),X2(end,2)+dose(2), X2(end,3:4)]);

    T = [T1(1:end);T2(2:end);T3(2:end)];
    X = [X1(1:end,1:3);X2(2:end,1:3);X3(2:end,1:3)];

    numsim = {X,T};

end
