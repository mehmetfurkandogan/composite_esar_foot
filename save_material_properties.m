%% Defining the material properties
clear;clc;close all;
t  =  0.131e-3; % m
rho = 1570;     % kg/m^3
E_1  =  177e9; % Pa          
E_2  =  12e9; % Pa          
E_3  =  12e9; % Pa          
G_12  = 5e9; % Pa          
G_31  = 5e9; % Pa          
G_23  = 4.4444e9; % Pa          
nu_12  =  0.3;
nu_13  =  0.3;
nu_23  =  0.35;

sigma_1_T_ult = 3310e6; % Pa
sigma_1_C_ult = 1793e6;  % Pa
sigma_2_T_ult = 96e6;   % Pa
sigma_2_C_ult = 240e6;  % Pa
tau_12_ult = 128e6;      % Pa

save('Materials/Ply_HexPly_8552_UD_IM10.mat')