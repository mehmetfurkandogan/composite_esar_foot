%% CLT Calculations for an ESR Foot
% Mehmet Furkan Doğan
% 30.11.2023
clc;clear;close all;
%% Defining the material properties
t  =  0.142e-3; % m
% Environment : T22C-m0% ; T = 22 °C ; m = 0 w%
% Engineering constants
E_1  =  156.5e9; % Pa          
E_2  =   8.83e9; % Pa          
E_3  =   8.83e9; % Pa          
G_12  =      4.3e9; % Pa          
G_31  =      4.3e9; % Pa          
G_23  =  3.39615e9; % Pa          
nu_12  =  0.3;
nu_13  =  0.3;
nu_23  =  0.3;

%% Calculation of compliance and stiffness matrices
% Compliance matrix for unidirectional lamina
S11 = 1/E_1;        % 1/Pa
S12 = -nu_12/E_1;   % 1/Pa
S13 = -nu_13/E_1;   % 1/Pa
S22 = 1/E_2;        % 1/Pa
S23 = -nu_23/E_2;   % 1/Pa
S33 = 1/E_3;        % 1/Pa
S44 = 1/G_23;       % 1/Pa
S55 = 1/G_31;       % 1/Pa
S66 = 1/G_12;       % 1/Pa
% Complete compliance matrix
S = [S11 S12 S13 0   0   0;
     S12 S22 S23 0   0   0;
     S13 S23 S33 0   0   0;
     0   0   0   S44 0   0;
     0   0   0   0   S55 0;
     0   0   0   0   0   S66];
% Stiffness matrix
C = inv(S);
% Reduced compliance matrix
S = [S11 S12 0;
     S12 S22 0;
     0   0   S66];  % 1/Pa
% Reduced stiffness matrix
Q = inv(S);         % Pa
%% Laminate Properties
% [45/-45/45/-45/0/90/0/90/-45/0/45/90]_s
theta = [45 -45 45 -45 0 90 0 90 -45 0 45 90];     % degree
theta = [theta flip(theta)];                       % degree (Symmetric)
n = size(theta,2);  % number of plies
H = n*t;        % m % Total width of the lamimate
for i = 0:n
    h(i+1) = -H/2 + i*t; % m
end

%% Angle transformation
% Reuter matrix
R = [1 0 0;
     0 1 0;
     0 0 2];    % -

for i = 1:n
    s = sin(theta(i)*pi/180);
    c = cos(theta(i)*pi/180);
    % Transformation matrix
    T = [c^2    s^2     2*s*c;
         s^2    c^2     -2*s*c;
         -s*c   s*c     c^2-s^2];  % -
    Qbar(:,:,i) = inv(T)*Q*R*T*inv(R);  % Pa
    Sbar(:,:,i) = inv(Qbar(:,:,i));     % 1/Pa
end

%% Global laminate stiffness matrix
A = zeros(3,3); % Pa * m    % Extensional stiffnes matrix
B = zeros(3,3); % Pa * m^2  % Coupling stiffness matrix
D = zeros(3,3); % Pa * m^3  % Bending stiffness matrix

for i = 1:3
    for j = 1:3
        for k = 1:n
            A(i,j) = A(i,j) + Qbar(i,j,k) * (h(k+1)-h(k));
            B(i,j) = B(i,j) + (1/2) * Qbar(i,j,k) * (h(k+1)^2-h(k)^2);
            D(i,j) = D(i,j) + (1/3) * Qbar(i,j,k) * (h(k+1)^3-h(k)^3);
        end
    end
end

ABBD = [A B;
        B D];
%% Loading
load('gait_forces.mat')
% Toe contact
[M1,I1] = max(F_foot_ground);
F_y = F_foot_ground_xp(I);
F_z = F_foot_ground_yp(I);
y_F = CoP_xp(I);
z_F = CoP_yp(I);
% Heel contact
[M2,I2] = max(F_foot_ground(1:46));