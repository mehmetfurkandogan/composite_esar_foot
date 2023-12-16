%% Defining the material properties
clear;clc;close all;
t  =  0.142e-3; % m
rho = 1583;     % kg/m^3
E_1  =  156.5e9; % Pa          
E_2  =  8.83e9; % Pa          
E_3  =    8.83e9; % Pa          
G_12  =      4.34e9; % Pa          
G_31  =      4.3e9; % Pa          
G_23  =      3.39615e9; % Pa          
nu_12  =  0.3;
nu_13  =  0.3;
nu_23  =  0.3;

sigma_1_T_ult = 2468e6; % Pa
sigma_1_C_ult = 1482e6;  % Pa
sigma_2_T_ult = 38e6;   % Pa
sigma_2_C_ult = 176.6e6;  % Pa
tau_12_ult = 128e6;      % Pa

save('Materials/Cycom 381 IM7 UD.mat')
%%
clear;clc;close all;
t_core  =  5e-3; % m
rho_core = 1495;     % kg/m^3
E_core  =  60.4e9; % Pa                   
nu_core  =  0.35;
G_core = E_core/(2*(1+nu_core));

save('Materials/Rohacell.mat')