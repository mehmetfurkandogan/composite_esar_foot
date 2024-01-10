function SR_inv = CLT(input)

t_core = input(end);
theta = input(1:end-1)*45;
stack = ceil(length(theta)/3);

%% Defining the material properties
core = true;

load('Materials/Cycom 381 IM7 UD.mat')

if t_core ~= 0
    load("Materials\Rohacell.mat")
    t_core = t_core*1e-3;
end
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
if t_core ~= 0
    theta_up = [theta(1:2*stack) 0 flip(theta(1:2*stack))];
    theta_down = [theta(2*stack+1:end) 0 flip(theta(2*stack+1:end))];  % degree (Symmetric)
    theta = [theta_up theta_down];
    n = size(theta,2);  % number of plies
    H = (n-2)*t+2*t_core;        % m % Total width of the lamimate
    h = zeros(1,n);
    h(1) = -H/2;
    for i = 1:length(theta)
        if i == 2*stack+1 || i == length(theta) + stack + 1
            h(i+1) = h(i) + t_core;
        else
            h(i+1) = h(i) + t; % m
        end
    end
else
    theta_up = [theta(1:2*stack) flip(theta(1:2*stack))];
    theta_down = [theta(2*stack+1:end) flip(theta(2*stack+1:end))]; % degree (Symmetric)
    theta = [theta_up theta_down];
    n = size(theta,2);  % number of plies
    H = n*t;        % m % Total width of the lamimate
    for i = 0:length(theta)
        h(i+1) = -H/2 + i*t; % m
    end
end

%% Angle transformation
% Reuter matrix
R = [1 0 0;
     0 1 0;
     0 0 2];    % -

for i = 1:n
    if i == 1+ (n-1)/2 && core == true
        Qbar(:,:,i) = [1/E_core -nu_core/E_core 0;
                       -nu_core/E_core 1/E_core 0;
                       0 0 1/G_core];
        Sbar(:,:,i) = inv(Qbar(:,:,i));     % 1/Pa
        continue
    end
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
% Foot Dimensions
shoe_size = 42;     % eu
L_data = 230e-3;    % m % Total Lenght of the foot
% L_model = 206e-3;  % m for proted design
L_model = ((shoe_size - 2 ) * 20 / 3)*1e-3;
b_rear = 0.311 * L_model;  % Lenght of the front part of the foot
b_front = L_model - b_rear;    % Lenght of the rear part of the foot
a = 0.31 * L_model; % Width of the foot
delta = 7.5;    % deg % the angle between foot axis and the walking direction
Iy = (1/12) * a*H^3; % m^4
Iz = (1/12) * H*a^3; % m^4
% J = Iy + Iz; % m^4
J = 1/16*a*H^3*(16/3-3.36*H/a*(1-H^4/(12*a^4)));

% Curved Beam
R = 100/208*L_model;
rn = H/log((R + H)/R);
e = R + H/2 - rn;
%% 
load('gait_forces.mat')

mass_data = 56.7;
design_mass = 130;

number_of_time_steps = length(spi);
nots = number_of_time_steps;
% Toe contact
Fz = -F_foot_ground_yp(spi)*design_mass/mass_data;
Fy = F_foot_ground_xp(spi) * sind(delta)*design_mass/mass_data;
Fx = F_foot_ground_xp(spi) * cosd(delta)*design_mass/mass_data;
b = CoP_xp(spi)*1e-3 * L_model / L_data - b_rear;
% b = abs(b);
%% OUTER FOOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loadings

Nx = Fx/(H*a)*H + Fy.*b*(a/2)*H/Iz + Fz.*b/(H*a*e)*(R*log((R+H/2+e)/(R-H/2+e)) - H);    % N/m
Ny = Fy./(H*b)*H;    % N/m
Nxy = Fy/(H*a)*H; %+ Fz*a^3/(8*J)*(H/a*sqrt(1+H^2/a^2) + 1/2*log(abs(H/a + sqrt(1+H^2/a^2))/abs(-H/a + sqrt(1+H^2/a^2))));  % N/m
N = [Nx Ny Nxy]';

Mx = Fz.*b/(2*H*a*e)*(2*H*R + 2*R*(R+e)*log((R-H/2+e)/(R+H/2+e)));   % N*m/m
My = zeros(size(b));      % N*m/m
Mxy = Fz*a^4/(16*J)*(-1/4*(H/a*sqrt(1+H^2/a^2) + 1/2*log(abs(H/a + sqrt(1+H^2/a^2))/abs(-H/a + sqrt(1+H^2/a^2)))) + 1/2*H/a*(H^2/a^2 + 1)^(3/2));                 % N*m/m
M = [Mx My Mxy]';

NM = [N;M];
%% Strains and curvatures
eps0kappa = ABBD\NM;

eps0 = eps0kappa(1:3,:);      % m/m
kappa = eps0kappa(4:6,:);     % 1/m

z = -H/2:100*1e-6:H/2;                     % z for whole laminate
for i = 1:nots
    eps(:,:,i) = repmat(eps0(:,i), 1, length(z)) + z .* kappa(:,i);
end
%% Stresses
sigma = zeros(size(eps));
sigma_loc = zeros(size(eps));
for i = 1:length(z)

    for j = 1:length(h)-1
        if z(i)>=h(j) && z(i)<h(j+1)
             ply = j;
             break
        end
    end
    if z(i)==h(end)
        ply = length(h)-1;
    end
    
    for j = 1:nots
        sigma(:,i,j) = Qbar(:,:,ply) * eps(:,i,j);
    end
    % Local stresses
    s = sin(theta(ply)*pi/180);
    c = cos(theta(ply)*pi/180);
    % Transformation matrix
    T = [c^2    s^2     2*s*c;
         s^2    c^2     -2*s*c;
         -s*c   s*c     c^2-s^2];  % -
    for j = 1:nots
        sigma_loc(:,i,j) = T * sigma(:,i,j);
    end
    
end

%% Strength Ratio
SR = zeros(length(z),nots);
% Tsai-Wu Criterion
H1 = 1/sigma_1_T_ult - 1/sigma_1_C_ult;     % 1/Pa
H11 = 1/(sigma_1_T_ult*sigma_1_C_ult);      % 1/Pa^2
H2 = 1/sigma_2_T_ult - 1/sigma_2_C_ult;     % 
H22 = 1/(sigma_2_T_ult*sigma_2_C_ult);
H6 = 0;
H66 = 1/tau_12_ult^2;
H12 = -1/2 * sqrt(1/(sigma_1_T_ult*sigma_1_C_ult*sigma_2_T_ult*sigma_2_C_ult)); % Mises-Hencky Criterion
for i=1:length(z)
    for j = 1:nots
        p = [H11*sigma_loc(1,i,j)^2+H22*sigma_loc(2,i,j)^2+H66*sigma_loc(3,i,j)^2+...
             H12*sigma_loc(1,i,j)*sigma_loc(2,i,j)...
             H1*sigma_loc(1,i,j)+H2*sigma_loc(2,i,j)+H6*sigma_loc(3,i,j) -1];
        SR(i,j) = max(roots(p));
    end
end

SR_tw = SR;

%Tsai-Hill Criterion
SR = zeros(length(z),nots);
for i=1:length(z)
    for j = 1:nots
        FI = sigma_loc(1,i,j)^2/sigma_1_T_ult^2 - sigma_loc(1,i,j)*sigma_loc(2,i,j)/sigma_1_T_ult^2 + sigma_loc(2,i,j)^2/sigma_2_T_ult^2 + sigma_loc(3,i,j)^2/tau_12_ult^2;
        SR(i,j) = 1/sqrt(FI);
    end
end

SR_th = SR;

%Max Stress Criterion
FI_1 = reshape(max(sigma_loc(1,:,:),0)/sigma_1_T_ult,[length(z),nots]);
FI_2 = reshape(max(-sigma_loc(1,:,:),0)/sigma_1_C_ult,[length(z),nots]);
FI_3 = reshape(max(sigma_loc(2,:,:),0)/sigma_2_T_ult,[length(z),nots]);
FI_4 = reshape(max(-sigma_loc(2,:,:),0)/sigma_2_C_ult,[length(z),nots]);
FI_5 = reshape(abs(sigma_loc(3,:,:))/tau_12_ult,[length(z),nots]);

SR_ms = 1./max(cat(3,FI_1, FI_2, FI_3, FI_4, FI_5),[],3);
SR_out = min(cat(3,SR_ms, SR_th, SR_tw),[],3);
%% INNER FOOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loadings

Nx = Fx/(H*a)*H + Fy.*b*(a/2)*H/Iz + Fz.*b/(H*a*e)*(R*log((R+H/2+e)/(R-H/2+e)) - H);    % N/m
Ny = -Fy./(H*b)*H;    % N/m
Nxy = -Fy/(H*a)*H; %+ Fz*a^3/(8*J)*(H/a*sqrt(1+H^2/a^2) + 1/2*log(abs(H/a + sqrt(1+H^2/a^2))/abs(-H/a + sqrt(1+H^2/a^2))));  % N/m
N = [Nx Ny Nxy]';

Mx = Fz.*b/(2*H*a*e)*(2*H*R + 2*R*(R+e)*log((R-H/2+e)/(R+H/2+e)));   % N*m/m
My = zeros(size(b));      % N*m/m
Mxy = -Fz*a^4/(16*J)*(-1/4*(H/a*sqrt(1+H^2/a^2) + 1/2*log(abs(H/a + sqrt(1+H^2/a^2))/abs(-H/a + sqrt(1+H^2/a^2)))) + 1/2*H/a*(H^2/a^2 + 1)^(3/2));                 % N*m/m
M = [Mx My Mxy]';

NM = [N;M];
%% Strains and curvatures
eps0kappa = ABBD\NM;

eps0 = eps0kappa(1:3,:);      % m/m
kappa = eps0kappa(4:6,:);     % 1/m

z = -H/2:100*1e-6:H/2;                     % z for whole laminate
for i = 1:nots
    eps(:,:,i) = repmat(eps0(:,i), 1, length(z)) + z .* kappa(:,i);
end
%% Stresses
sigma = zeros(size(eps));
sigma_loc = zeros(size(eps));
for i = 1:length(z)

    for j = 1:length(h)-1
        if z(i)>=h(j) && z(i)<h(j+1)
             ply = j;
             break
        end
    end
    if z(i)==h(end)
        ply = length(h)-1;
    end
    
    for j = 1:nots
        sigma(:,i,j) = Qbar(:,:,ply) * eps(:,i,j);
    end
    % Local stresses
    s = sin(theta(ply)*pi/180);
    c = cos(theta(ply)*pi/180);
    % Transformation matrix
    T = [c^2    s^2     2*s*c;
         s^2    c^2     -2*s*c;
         -s*c   s*c     c^2-s^2];  % -
    for j = 1:nots
        sigma_loc(:,i,j) = T * sigma(:,i,j);
    end
    
end

%% Strength Ratio
SR = zeros(length(z),nots);
% Tsai-Wu Criterion
H1 = 1/sigma_1_T_ult - 1/sigma_1_C_ult;     % 1/Pa
H11 = 1/(sigma_1_T_ult*sigma_1_C_ult);      % 1/Pa^2
H2 = 1/sigma_2_T_ult - 1/sigma_2_C_ult;     % 
H22 = 1/(sigma_2_T_ult*sigma_2_C_ult);
H6 = 0;
H66 = 1/tau_12_ult^2;
H12 = -1/2 * sqrt(1/(sigma_1_T_ult*sigma_1_C_ult*sigma_2_T_ult*sigma_2_C_ult)); % Mises-Hencky Criterion
for i=1:length(z)
    for j = 1:nots
        p = [H11*sigma_loc(1,i,j)^2+H22*sigma_loc(2,i,j)^2+H66*sigma_loc(3,i,j)^2+...
             H12*sigma_loc(1,i,j)*sigma_loc(2,i,j)...
             H1*sigma_loc(1,i,j)+H2*sigma_loc(2,i,j)+H6*sigma_loc(3,i,j) -1];
        SR(i,j) = max(roots(p));
    end
end

SR_tw = SR;

%Tsai-Hill Criterion
SR = zeros(length(z),nots);
for i=1:length(z)
    for j = 1:nots
        if sigma_loc(1,i,j) < 0 X1 = sigma_1_C_ult; else X1 = sigma_1_T_ult; end
        if sigma_loc(2,i,j) < 0 X2 = sigma_2_C_ult; else X2 = sigma_2_T_ult; end
        FI = sigma_loc(1,i,j)^2/X1^2 - sigma_loc(1,i,j)*sigma_loc(2,i,j)/X1^2 + sigma_loc(2,i,j)^2/X2^2 + sigma_loc(3,i,j)^2/tau_12_ult^2;
        SR(i,j) = 1/sqrt(FI);
    end
end

SR_th = SR;

%Max Stress Criterion
FI_1 = reshape(max(sigma_loc(1,:,:),0)/sigma_1_T_ult,[length(z),nots]);
FI_2 = reshape(max(-sigma_loc(1,:,:),0)/sigma_1_C_ult,[length(z),nots]);
FI_3 = reshape(max(sigma_loc(2,:,:),0)/sigma_2_T_ult,[length(z),nots]);
FI_4 = reshape(max(-sigma_loc(2,:,:),0)/sigma_2_C_ult,[length(z),nots]);
FI_5 = reshape(abs(sigma_loc(3,:,:))/tau_12_ult,[length(z),nots]);

SR_ms = 1./max(cat(3,FI_1, FI_2, FI_3, FI_4, FI_5),[],3);
SR_in = min(cat(3,SR_ms, SR_th, SR_tw),[],3);
%% Output
SR = min(min(SR_in,SR_out),[],'all');
SR_inv = 1/SR;
end