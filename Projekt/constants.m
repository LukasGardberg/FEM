% Loads constants for FEM program

% SUBDOMAINS
%
% 1: Silver-epoxy
% 2 & 3: Copper
% 4: Silicon


% Material densities
ro_epoxy = 2500;
ro_silicon = 2530;
ro_copper = 8930;

ro_const= [ro_epoxy, ro_copper, ro_copper, ro_silicon];


% Youngs modulus, E

e_epoxy = 7;
e_silicon = 165;
e_copper = 128;

e_const = [e_epoxy, e_copper, e_copper, e_silicon]*10^9;


% Poissons ratio, ny

ny_epoxy = 0.3;
ny_silicon = 0.22;
ny_copper = 0.36;

ny_const = [ny_epoxy, ny_copper, ny_copper, ny_silicon];


% Expansion coefficient, alpha

alpha_epoxy = 40*10^(-6);
alpha_silicon = 2.6*10^(-6);
alpha_copper = 17.6 * 10^(-6);

alpha_const = [alpha_epoxy, alpha_copper, alpha_copper, alpha_silicon];


% Conductivity constants

k_epoxy = 5;
k_silicon = 149;
k_copper = 385;

k_const = [k_epoxy, k_copper, k_copper, k_silicon];


% Specific heat

c_epoxy = 1000;
c_silicon = 703;
c_copper = 386;

c_const = [c_epoxy, c_copper, c_copper, c_silicon];


% Thickness: 10mm
thickness = 10 * 10^(-3); % in meters

% Q = 5 * 10^7; % W / m^3, r??tt enhet??

T_0 = 30;
T_inf = 18;