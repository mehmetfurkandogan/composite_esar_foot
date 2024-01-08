%% Optimization of layer orientations
% Optimizations are done using genetic algorithm
% 16.12.2023
clc;clear;close all;
%%
laminaMin = 7*4;
laminaMax = 10*4;

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

for iter = 1:(laminaMax-laminaMin)/4 + 1
    laminaCount = 4*(iter-1) + laminaMin;
    lb = ones(1, laminaCount);
    ub = zeros(1, laminaCount);
    for i = 1:laminaCount
        ub(i) = laminaCount-i+1;
    end
    tic
    [Tbest{iter}, SRbest(iter), ~] = ga(@CLT, ...
        laminaCount, [], [], [], [], lb, ub, [], 1:laminaCount, opts);
    toc
    options = [zeros(1,laminaCount/4) 45*ones(1,laminaCount/4) -45*ones(1,laminaCount/4) 90*ones(1,laminaCount/4)];

    theta = [];
    for i = 1:laminaCount
        theta = [theta options(Tbest{iter}(i))];
        options(Tbest{iter}(i)) = [];
    end
    fprintf('Optimum layer orientation for %d laminates:\t',laminaCount);
    disp(theta);
    fprintf('Maximum strength ratio for %d laminates:\t%.2f\n',laminaCount,1/SRbest(iter));
end