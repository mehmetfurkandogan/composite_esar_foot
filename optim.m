clc;
clear;

laminaMin = 5;
laminaMax = 6;

Tbest = cell(1, laminaMax - laminaMin + 1);
SRbest = zeros(1, laminaMax - laminaMin + 1);
exitflag = zeros(1, laminaMax - laminaMin + 1);

opts = optimoptions(@ga, ...
                    'PopulationSize', 150, ...
                    'MaxGenerations', 200, ...
                    'EliteCount', 10, ...
                    'FunctionTolerance', 1e-8, ...
                    'PlotFcn', @gaplotbestf);

for iter = 1:laminaMax-laminaMin + 1
    laminaCount = iter + laminaMin - 1;
    lb = -ones(1, laminaCount);
    ub = 2*ones(1, laminaCount);
    [Tbest{iter}, SRbest(iter), ~] = ga(@CLT, ...
        laminaCount, [], [], [], [], lb, ub, [], 1:laminaCount, opts);
end