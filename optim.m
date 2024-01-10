%% Optimization of layer orientations
% Optimizations are done using genetic algorithm
% 16.12.2023
clc;clear;close all;
%%
laminaMin = 3*12;
laminaMax = 3*15;

Tbest = cell(1, laminaMax - laminaMin + 1);
SRbest = zeros(1, laminaMax - laminaMin + 1);
exitflag = zeros(1, laminaMax - laminaMin + 1);

opts = optimoptions(@ga, ...
                    'PopulationSize', 10, ...
                    'MaxGenerations', 200, ...
                    'EliteCount', 2, ...
                    'FunctionTolerance', 1e-6, ...
                    'PlotFcn', @gaplotbestf,...
                    'Display','iter');

for iter = 1:(laminaMax-laminaMin)/3 + 1
    laminaCount = 3*(iter-1) + laminaMin;
    lb = ones(1, laminaCount);
    ub = zeros(1, laminaCount);

    for i = 1:2*laminaCount/3
        ub(i) = (2*laminaCount/3)-i+1;
    end

    for i = 2*laminaCount/3+1:laminaCount
        ub(i) = (laminaCount/3)-i+2*laminaCount/3+1;
    end

    tic
    [Tbest{iter}, SRbest(iter), ~] = ga(@CLT, ...
        laminaCount, [], [], [], [], lb, ub, [], 1:laminaCount, opts);
    toc

    % Upper Kneel
    stack = laminaCount/3;
    limit = ceil(stack*2/8);
    remaining = 2*stack - 3*limit;
    
    %UD
    options = [zeros(1,remaining) 45*ones(1,limit) -45*ones(1,limit) 90*ones(1,limit)];
    
    %Woven
    % options = [zeros(1,length(layers)/2) 45*ones(1,length(layers)/2)];
    
    theta_up = [];
    for i = 1:2*stack
        theta_up = [theta_up options(Tbest{iter}(i))];
        options(Tbest{iter}(i)) = [];
    end

    stack = laminaCount/3;
    limit = ceil(stack/8);
    remaining = stack - 3*limit;
    
    %UD
    options = [zeros(1,remaining) 45*ones(1,limit) -45*ones(1,limit) 90*ones(1,limit)];
    
    %Woven
    % options = [zeros(1,length(layers)/2) 45*ones(1,length(layers)/2)];
    
    theta_low = [];
    for i = 2*stack+1:laminaCount
        theta_low = [theta_low options(Tbest{iter}(i))];
        options(Tbest{iter}(i)) = [];
    end

    fprintf('UPPER KNEEL: Optimum layer orientation for %d laminates:\t',laminaCount);
    disp(theta_up);
    fprintf('LOWER KNEEL: Optimum layer orientation for %d laminates:\t',laminaCount);
    disp(theta_low);
    fprintf('Maximum strength ratio for %d laminates:\t%.2f\n',laminaCount,1/SRbest(iter));
end