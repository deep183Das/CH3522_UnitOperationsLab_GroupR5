%% CH3522 Experiment 6: Shell-Tube Heat Exchanger
%  Author: Deepanjhan Das [CH22B020]
clear; clc; close all;
format long;

%% Experimental Data 

% hot water V flow rate (LPH) -> (m^3/s)
Vhot = [ones(1,3)*50, ones(1,3)*125, ones(1,3)*175] * (0.001 / 3600);
% cold water V flow rate (LPH) -> (m^3/s)
Vcool = [70, 155, 195];
Vcold = [Vcool, Vcool, Vcool] * (0.001 / 3600);

% hot water inlet temperature (K)
T1 = [61.1, 60.2, 59.7, 61.5, 60.2, 61.0, 60.5, 59.5, 60.6] + 273.15; 
% hot water outlet temperature (K)
T2 = [42.6, 47.3, 39.0, 46.1, 44.2, 45.3, 50.1, 47.6, 47.5] + 273.15;

% cold water inlet temperature (K)
T3 = [29.5, 29.9, 29.7, 29.7, 30.1, 30.1, 29.7, 30.0, 30.0] + 273.15;
% cold water outlet temperature (K)
T4 = [39.9, 43.7, 35.9, 43.2, 40.1, 40.0, 45.9, 42.4, 41.5] + 273.15;


%% Provided properties

% hot water density (kg/m^3)
rho_hot = 982.60;
% cold water density (kg/m^3)
rho_cold = 995.71;

% hot water specific heat capacity (J/kg.K) @ 60 deg C
Cp_hot = 4185;
% cold water specific heat capacity (J/kg.K) @ 30 deg C
Cp_cold = 4178;

% Number of tubes
N = 24;
% inner diameter (m)
Di = 0.013;
% outer diameter (m)
Do = 0.016;
% length of each tube (m)
L = 0.5;


%% Calculations I (Includes 3 Objectives).
% -> Calculation of the amount of the heat transferred.

% Function for amount of heat gained (W or J/s)
Q = @(rho, V, Cp, Ti, To) (rho .* V .* Cp .* (To - Ti));

% heat gained by cold water per s (W)
Qcold = Q(rho_cold, Vcold, Cp_cold, T3, T4);

% heat lost by hot water per s (W)
Qhot = abs(Q(rho_hot, Vhot, Cp_hot, T1, T2));

% average of the heat lost & gained per s (W)
Qavg = (Qcold + Qhot) ./ 2;
fprintf('Q Calcualtions Done!\n');

% -> Del_T_LMTD calculation.

% Function for "log-mean-temperature difference)
Del_T_LMTD = @(t1, t2, t3, t4) (((t1 - t3) - (t2 - t4)) ./ log((t1 - t3)./(t2 - t4)));

% LMTD (K)
LMTD = Del_T_LMTD(T1, T2, T3, T4);
fprintf('LMTD Calculation Done!\n');

% -> Calculation of Overall HTC (W/m^2.K).

% inner area (m^2)
Ai = pi * Di * L * N;
% outer area (m^2)
Ao = pi * Do * L * N;

% In a real scenario Qcold & Qhot won't be equal due the inefficiency in
% the system, i.e. the heat exchanger. So we use the averaged value in the
% calculation of overall HTC both at inner and outer area.

% Function for HTC (inner & outer) (W/m^2.K)
U = @(A) (Qavg ./ (A .* LMTD));

% inner surface Overall HTC
Ui = U(Ai);
% outer surface Overall HTC
Uo = U(Ao);
fprintf("Calculation for U Done!\n");


%% Calculations II (Calculation of Theoretical Overall HTC)

% inner cross-sectional area (m^2)
ai = pi * (Di/2)^2;
% outer cross-sectional area (m^2)
ao = pi * (Do/2)^2;

% velocity of hot water (m/s)
vhot = Vhot / ai;
% velocity of cold water (m/s)
vcold = Vcold / ao;

% dynamic viscosity of hot water (Pa.s)
mu_hot = 0.0004656;
% dynamic viscosity of cold water (Pa.s)
mu_cold = 0.000797;

% thermal conductivity of hot water (W/m.K)
Khot = 0.65091;
% thermal conductivity of cold water (W/m.K)
Kcold = 0.61450;

% Functions for the dimensionless numbers.
% Reynold's number
Re = @(rho, v, D, mu) ((rho .* v .* D) ./ mu);

% Prandtl number 
Pr = @(mu, Cp, K) ((mu * Cp) / K);

% Heat Transfer Coefficients 
h = @(NuN, K, D) ((NuN .* K) ./ D);

% inner surface theoretical results
Rei = Re(rho_hot, vhot, Di, mu_hot);
Pri = Pr(mu_hot, Cp_hot, Khot);
Nui = 0.23 .* (Rei.^0.8) .* (Pri.^0.3); % Nusselt number for cooling
hi = h(Nui, Khot, Di);

% outer surface theoretical results
Reo = Re(rho_cold, vcold, Do, mu_cold);
Pro = Pr(mu_cold, Cp_cold, Kcold);
Nuo = 0.23 .* (Reo.^0.8) .* (Pro.^0.4); % Nusselt number for heating
ho = h(Nuo, Kcold, Do);

% Overall HTC (W/m^2.K)
% Considering that Carbon Steel is the most standard material to be used as
% the pipe material, we need its thermal conductivity value. 
K = 54; % in (W/m.K)
U_theoretical = ((1./ho) + Do./(hi.*Di) + Do.*((log(Do./Di))./(2*K))).^(-1);
fprintf("Theoretical Results' Calculations Done!\n");