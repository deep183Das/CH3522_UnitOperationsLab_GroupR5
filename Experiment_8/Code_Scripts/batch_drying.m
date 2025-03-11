%% CH3522 Experiment 8: Batch Drying
%  Author: Deepanjhan Das [CH22B020]
clear; clc; close all;
format long;

%% Data using the Experimental Tabulations.

% The Humid ratio is found using the online-calculator from the following,
% https://www.flycarpet.net/en/psyonline

% flow rates (LPM) -> m^3 / s
Q1 = 2/(60 * 1000);
Q2 = 6/(60 * 1000);
Q3 = 9/(60 * 1000);

% I am not mentioning the RH and Wet Bulb temperature again. Rather using
% the link provided above I am directly mentioning the Humid Ratio (HR) &
% Dry Bulb T, which will be used in further calculations. 
%
% HR in (g H2O / Kg Dry Air)
% Tdb in K
% Psat in Pa
HR1 = [17.796, 17.397, 17.059, 16.856, 16.769, 16.778, 16.682, 16.821, ...
       17.019, 17.332, 17.357, 17.357, 17.349, 17.357, 17.357, 17.506, 17.506, ...
       17.697, 17.832, 17.673, 17.539, 17.556, 17.564, 17.564];

Tdb1 = [31.760, 31.520, 31.154, 30.862, 30.685, 30.665, 30.510, 30.564, ...
        30.477, 30.506, 30.447, 30.447, 30.466, 30.447, 30.447, 30.482, 30.482, ...
        30.419, 30.493, 30.477, 30.403, 30.365, 30.345, 30.345] + 273.15;

Psat1 = [4.696, 4.633, 4.537, 4.463, 4.418, 4.413, 4.374, 4.387, ...
         4.366, 4.373 , 4.358, 4.358, 4.363, 4.358, 4.358, 4.367, 4.367, ...
         4.351, 4.37, 4.366, 4.347, 4.338, 4.333, 4.333] .* 1000;

HR2 = [19.397, 19.622, 19.923, 19.960, 19.825, 19.810, 19.795, 19.675, 19.660, 19.660];
Tdb2 = [30.389, 30.265, 30.376, 30.288, 30.196, 30.232, 30.267, 30.140, 30.176, 30.176] + 273.15;
Psat2 = [4.344, 4.313, 4.34, 4.319, 4.296, 4.305, 4.313, 4.282, 4.291, 4.291] .* 1000;

HR3 = [21.432, 21.240, 20.893, 20.542, 20.157, 19.488, 19.256, 18.776, 18.615, 18.303, 18.287, 18.287];
Tdb3 = [30.172, 30.199, 30.170, 30.160, 30.238, 30.173, 30.315, 30.240, 30.222, 30.166, 30.204, 30.204] + 273.15; 
Psat3 = [4.29, 4.297, 4.289, 4.287, 4.306, 4.29, 4.325, 4.307, 4.302, 4.288, 4.298, 4.298] .* 1000;


%% Calculations for Air at inlet.

% (inlet condition same for all flow rates) RH = 52.8% & WB T = 22.9 deg C
HR_in = 14.518;
Tdb_in = 30.462 + 273.15;

MW_h2o = 18; % g/mol
MW_air = 28.96;

% mole ratio (H2O / D.A.)
moleRatio_in = (HR_in .* MW_air) ./ (1000 .* MW_h2o);
% mole fraction (H2O / (H2O + D.A))
xin = moleRatio_in ./ (1 + moleRatio_in);

% Atmospheric Pressure (Pa)
P_atm = 101325;

% Vapor Pressure at inlet (Pa)
VP_in = xin .* P_atm;

% Assuming ideal gas law to be valid, molar flow rate of
% water vapor at inlet is as follows, (mol / s).
R = 8.314;
molar_flow_rate_in_1 = (VP_in * Q1) ./ (R .* Tdb_in);
molar_flow_rate_in_2 = (VP_in * Q2) ./ (R .* Tdb_in);
molar_flow_rate_in_3 = (VP_in * Q3) ./ (R .* Tdb_in);


%% Calculations for Air at outlet.

% mole ratio (H2O / D.A.)
moleRatio1 = (HR1 .* MW_air) ./ (1000 .* MW_h2o);
moleRatio2 = (HR2 .* MW_air) ./ (1000 .* MW_h2o);
moleRatio3 = (HR3 .* MW_air) ./ (1000 .* MW_h2o);

% mole fraction (H2O / (H2O + D.A))
x1 = moleRatio1 ./ (1 + moleRatio1);
x2 = moleRatio2 ./ (1 + moleRatio2);
x3 = moleRatio3 ./ (1 + moleRatio3);

% Vapor Pressure at outlet (Pa)
VP1_out = x1 .* P_atm;
VP2_out = x2 .* P_atm;
VP3_out = x3 .* P_atm;

% molar flow rate (mol/s)
molar_flow_rate_out_1 = (VP1_out .* Q1) ./ (R .* Tdb1);
molar_flow_rate_out_2 = (VP2_out .* Q2) ./ (R .* Tdb2);
molar_flow_rate_out_3 = (VP3_out .* Q3) ./ (R .* Tdb3);


%% Other calculations.

% rate of evaporation (mol/s)
ev_rate_1 = molar_flow_rate_out_1 - molar_flow_rate_in_1;
ev_rate_2 = molar_flow_rate_out_2 - molar_flow_rate_in_2;
ev_rate_3 = molar_flow_rate_out_3 - molar_flow_rate_in_3;

% area of evaporation (m^2)
ev_area = 0.000625;

% molar flux (mol/m^2.s)
N1 = (ev_rate_1 ./ ev_area);
N2 = (ev_rate_2 ./ ev_area);
N3 = (ev_rate_3 ./ ev_area);

% N = Ky (Ys - Y)
%
% Y : average mole fraction
Y1 = (x1 + xin) ./ 2;
Y2 = (x2 + xin) ./ 2;
Y3 = (x3 + xin) ./ 2;

% Ys : mole fraction at saturation vapor pressure
% Saturation vapor pressure will come from the Psychrometric chart for Tdb.
Ys1 = Psat1 ./ P_atm;
Ys2 = Psat2 ./ P_atm;
Ys3 = Psat3 ./ P_atm;

% Ky (mol/m^2.s) using the above equation
Ky1 = N1 ./ (Ys1 - Y1);
Ky2 = N2 ./ (Ys2 - Y2);
Ky3 = N3 ./ (Ys3 - Y3);

% Volumetric Flow rate (mol/m^3)
VF1 = Psat1 ./ (R .* Tdb1);
VF2 = Psat2 ./ (R .* Tdb2);
VF3 = Psat3 ./ (R .* Tdb3);

% Kc (m/s) = Ky / VF 
Kc1 = Ky1 ./ VF1;
Kc2 = Ky2 ./ VF2;
Kc3 = Ky3 ./ VF3;


%% Mass Transfer Coefficient (Theoretical Correlation)
l = 0.025; % (m) characteristic length
Ac = pi * (0.006 / 2)^2; % (m^2) cross-sectional area

% velocity of air (m/s)
v1 = Q1 / Ac;
v2 = Q2 / Ac;
v3 = Q3 / Ac;

% Reynold's number
rho_air = 1.293; % kg / m^3
mu_air = 1.598 * 10^(-5); % m^2/s

N_Re1 = (rho_air * v1 * l) / mu_air;
N_Re2 = (rho_air * v2 * l) / mu_air;
N_Re3 = (rho_air * v3 * l) / mu_air;

% Schmidt number
D_air = 2.26 * 10^(-5); % m^2/s

N_Sc = mu_air / (rho_air * D_air);

% Mass transfer coefficient (m/s)
k_m1 = (D_air/l) * 0.646 * (N_Re1^0.5) * (N_Sc^(1/3));
k_m2 = (D_air/l) * 0.646 * (N_Re2^0.5) * (N_Sc^(1/3));
k_m3 = (D_air/l) * 0.646 * (N_Re3^0.5) * (N_Sc^(1/3));

% Show in figure
figure();
hold on;
grid on;
plot(1:length(Kc1), Kc1, 'r-', LineWidth=1.1, Marker='o', MarkerFaceColor='blue');
plot(1:length(Kc1), ones(1,size(Kc1,2))*k_m1, 'b-', LineWidth=1.1);
legend('Experimental Value', "Theoretical Value");
xlabel("Time (minutes)");
ylabel("Mass Transfer Coefficient (m/s)");
title("MTC for 2LPM flow rate of dry air");
hold off;

figure();
hold on;
grid on;
plot(1:length(Kc2), Kc2, 'r-', LineWidth=1.1, Marker='o', MarkerFaceColor='blue');
plot(1:length(Kc2), ones(1,size(Kc2,2))*k_m2, 'b-', LineWidth=1.1);
legend('Experimental Value', "Theoretical Value");
xlabel("Time (minutes)");
ylabel("Mass Transfer Coefficient (m/s)");
title("MTC for 6LPM flow rate of dry air");
hold off;

figure();
hold on;
grid on;
plot(1:length(Kc3), Kc3, 'r-', LineWidth=1.1, Marker='o', MarkerFaceColor='blue');
plot(1:length(Kc3), ones(1,size(Kc3,2))*k_m3, 'b-', LineWidth=1.1);
legend('Experimental Value', "Theoretical Value");
xlabel("Time (minutes)");
ylabel("Mass Transfer Coefficient (m/s)");
title("MTC for 9LPM flow rate of dry air");
hold off;
