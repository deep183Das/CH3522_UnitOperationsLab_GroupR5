%% CH3522 Experiment 6: Fin-Tube Heat Exchanger
%  Author: Deepanjhan Das [CH22B020]
clear; close all;
format long;

%% Experimental Data

% pressure of the steam (Pa)
P = [1, 1, 1, 2, 2, 2, 3, 3, 3] .* 6894.76;

% air h1 (m) & h2 (m)
h1 = [38.3, 36.4, 37.6, 37.5, 38.5, 38.2, 38.0, 38.9, 39.4] ./ 100;
h2 = [27.0, 29.4, 28.0, 28.0, 27.0, 27.3, 27.4, 26.7, 26.0] ./ 100;

% steam inlet temperature (K)
T1 = [118, 116, 117, 130, 130, 130, 140, 140, 141] + 273.15;
% steam outlet temperature (K)
T2 = [115, 102, 100, 130, 130, 130, 135, 122, 140] + 273.15;

% air inlet temperature (K)
T3 = [29, 29, 29, 28, 28, 28, 28, 28, 28] + 273.15;
% air outlet temperature (K)
T4 = [96, 99, 101, 112, 113, 115, 116, 120, 122] + 273.15;


%% Given data

% inner diameter of pipe (m)
pipe_id = 0.038;
% outer diameter of pipe (m)
pipe_od = 0.048;
% diameter of the orifice (m)
orifice_d = 0.025;


%% Experimental Determinations

% -> Flow rate of air
beta = orifice_d / pipe_id;
Cd = 0.9;
C = Cd / (sqrt(1 - (beta^4)));
A2 = pi * (orifice_d/2)^2;
rho_air = 1.164; % air density at 30 deg C (kg/m^3)

del_h = h1 - h2;
del_P = rho_air .* 9.8 .* del_h;
qm = rho_air .* C .* A2 .* (sqrt((2.*(del_P)).*(rho_air))); % mass flow rate (kg/s)


% -> Heat Transfer Area
% calculated separately.


% -> LMTD
% steam out - air in
delT1 = T2 - T3;
% steam in - air out
delT2 = T1 - T4;
LMTD = (delT1 - delT2) ./ (log((delT1)./(delT2)));


% -> U experimental

% wrt air (J/kg.K)
Cp_air = 1008; 
% heat transfered per time (W)
Q = qm .* Cp_air .* (T4 - T3);
% from the specified calculation value of A (m^2) comes out to be,
A = 0.47933;

% overall HTC
U = Q ./ (A .* LMTD);


%% Theoretical Determinations

% -> heating of Air
K_air = 0.026;
mu_air = 3.178 * 10^-5;
v_air = qm ./ (rho_air * A2);

Re_air = (rho_air .* v_air .* pipe_od) ./ mu_air;
Pr_air = (mu_air * Cp_air) ./ K_air;
Nu_air = 0.023 .* (Re_air.^0.8) .* (Pr_air.^0.4);
h_air = (Nu_air .* K_air) ./ pipe_od;


% -> cooling of steam
K_steam = 2150;
mu_steam = 3 * 10^-5;
rho_steam = 0.6;
Cp_steam = 1996;
v_steam = (A2 .* v_air) ./ (pi .* (pipe_id/2).^2);

Re_steam = (rho_steam .* v_steam .* pipe_id) ./ mu_steam;
Pr_steam = (mu_steam * Cp_steam) ./ K_steam;
Nu_steam = 0.023 .* (Re_steam.^0.8) .* (Pr_steam.^0.3);
h_steam = (Nu_steam .* K_steam) ./ pipe_id;


% thermal conductivity of the metal
K = 24.35;
A_fins = 0.265329;

% -> Fin efficiency
L = 1.49;
m = 110.62;
eta = tanh(m*L) / (m*L);
eta_f = 1 - (A_fins/A)*(1-eta);

% -> U theoretical
U_theo = ((pipe_od ./ (eta_f.*h_steam.*pipe_id)) + (pipe_od .* (log(pipe_od/pipe_id))) ./ (K) + 1./(eta_f.*h_air)).^-1;