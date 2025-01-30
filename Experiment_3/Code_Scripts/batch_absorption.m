%% CH 3522 Lab | Batch Absorption
%  Author: Deepanjhan Das [CH22B020]
format long;
clear; clc; close all;

% conversion factors
kelvin_factor = 273.15;
bar_to_pascal = 1e5;
mL_to_m3 = 1e-6;

V = 135; % volume (mL) of the vessel
R = 8.314; % J/mol.K


%% Properties of CO2
omega = 0.228; % ecentric factor
Tc = 304.25;   % critical temperature (K)
Pc = 7.39 * 1e6; % critical pressure (Pa)
kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega^2;
a = (0.45724 * R^2 * Tc^2) / Pc;
b = (0.07780 * R * Tc) / Pc;

% Defining anonymous functions
Tr = @(T) (T ./ Tc);
alpha = @(T) (1 + kappa.*(1 - (Tr(T)).^0.5)).^2;

% Defining the ideal state-equation
PIdeal = @(C, T) (C .* R .* T);
CIdeal = @(P, T) (P ./ (R.*T));
% Defining the real state-equation (Peng Robinson)
PReal = @(C, T) (((R*T)./((1./C)-b)) - ((a.*alpha(T))./((1./C.^2) + 2*b./C - b)));

% Define the function for fsolve (to get CReal for given P & T with an initial value of C)
pressure_equation = @(C, P, T) ...
    ((R.*T./((1./C)-b)) - ((a.*alpha(T))./((1./C.^2) + 2*b./C - b))) - P;


%% [Part 1] EtOH & CO2

% temp (deg C)
T1 = 27.0 + kelvin_factor; % [26.9 - 27.1]
% rpm (280-320)
rpm1 = 300; 
% time (min)
t1 = 0:1:31;
% pressure (bar)
P1 = [3.15, 2.30, 1.90, 1.70, 1.55, 1.45, 1.40, 1.30, 1.30, 1.25, 1.20];
P1 = [P1, ones(1,9)*1.15, ones(1,12)*1.10];

% Non-Linear regression is performed (using OLS objective)
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt', 'Display','iter');
% [Ideal Data] ------------------------------------------------------------
C1Idealexp = CIdeal(P1*bar_to_pascal, T1);
y11 = C1Idealexp;
x10(1) = 50;
x10(2) = 40;
x10(3) = 0.02;

x1IdealEst = lsqnonlin(@(x1)ConcFunc(x1,y11,t1), x10, [],[], options);

% Deriving the results
C1Idealmodel = x1IdealEst(1) + x1IdealEst(2) .* exp(-x1IdealEst(3) .* t1);
P1Idealest = PIdeal(C1Idealmodel, T1);


% [Real Data] -------------------------------------------------------------
C0 = 50; % mol / m^3
C1Realexp = zeros(size(P1));
for i=1:length(P1)
    C1Realexp(i) = fsolve(@(C) pressure_equation(C, P1(i)*bar_to_pascal, T1), C0);
end
y12 = C1Realexp;

x1RealEst = lsqnonlin(@(x1)ConcFunc(x1,y12,t1), x10, [],[], options);

% Deriving the results
C1Realmodel = x1RealEst(1) + x1RealEst(2) .* exp(-x1RealEst(3) .* t1);
P1Realest = PReal(C1Realmodel, T1);


%% Visualize
toSave = 'none';
% Concentration values [Ideal model]
y_label = sprintf('Concentration (mol/m^3)');
title_name = sprintf('Concentration vs Time for EtOH - CO_2 mixture [Ideal State Model]');
legend_label = sprintf('Ideal State est');
% toSave = 'C_etohco2_ideal.png';
text_label = sprintf('C(mol/m^3)[t] = %.4f + (%.4f)*exp(-%.4f*t)', x1IdealEst(1), x1IdealEst(2), x1IdealEst(3));
visualize1(C1Idealexp, C1Idealmodel, t1, title_name, y_label, legend_label, text_label, 10, 98, toSave);

% Concentration values [Ideal model]
y_label = sprintf('Concentration (mol/m^3)');
title_name = sprintf('Concentration vs Time for EtOH - CO_2 mixture [Real State Model]');
legend_label = sprintf('Real State est');
% toSave = 'C_etohco2_real.png';
text_label = sprintf('C(mol/m^3)[t] = %.4f + (%.4f)*exp(-%.4f*t)', x1RealEst(1), x1RealEst(2), x1RealEst(3));
visualize1(C1Realexp, C1Realmodel, t1, title_name, y_label, legend_label, text_label, 10, 98, toSave);

% Pressure values
y_label = sprintf('P (bar)');
title_name = sprintf('Pressure vs Time for EtOH - CO_2 mixture [Both Models]');
% toSave = 'P_etohco2.png';
visualize2(P1, (P1Idealest./bar_to_pascal), (P1Realest./bar_to_pascal), t1, title_name, y_label, toSave);


%% [Part 2] H2O & CO2

% temp (deg C)
T2 = 26.7 + kelvin_factor; % [26.7 - 26.8]
% rpm (310-330)
rpm2 = 320; 
% time (min)
t2 = 0:0.5:3; 
t2 = [t2, 4:1:58];
% pressure (bar)
P2 = [3.10, 2.95, 2.95, 2.90, 2.85, 2.85, 2.85, 2.80, 2.75, 2.70, 2.70,...
      2.60, 2.60, 2.55];
P2 = [P2, ones(1,3)*2.50, 2.45, ones(1,3)*2.40, 2.35, ones(1,3)*2.30, 2.25,...
      ones(1,2)*2.20, ones(1,5)*2.15, ones(1,2)*2.10, ones(1,5)*2.05, ...
      ones(1,2)*2.0, ones(1,6)*1.95, ones(1,3)*1.90, ones(1,11)*1.85];

% Non-Linear regression is performed (using OLS objective)
% [Ideal Data] ------------------------------------------------------------
C2Idealexp = CIdeal(P2*bar_to_pascal, T2);
y21 = C2Idealexp;
x20(1) = 50;
x20(2) = 40;
x20(3) = 0.02;

x2IdealEst = lsqnonlin(@(x2)ConcFunc(x2,y21,t2), x20, [],[], options);

% Deriving the results
C2Idealmodel = x2IdealEst(1) + x2IdealEst(2) .* exp(-x2IdealEst(3) .* t2);
P2Idealest = PIdeal(C2Idealmodel, T2);


% [Real Data] -------------------------------------------------------------
C2Realexp = zeros(size(P2));
for i=1:length(P2)
    C2Realexp(i) = fsolve(@(C) pressure_equation(C, P2(i)*bar_to_pascal, T2), C0);
end
y22 = C2Realexp;

x2RealEst = lsqnonlin(@(x2)ConcFunc(x2,y22,t2), x20, [],[], options);

% Deriving the results
C2Realmodel = x2RealEst(1) + x2RealEst(2) .* exp(-x2RealEst(3) .* t2);
P2Realest = PReal(C2Realmodel, T2);


%% Visualize the estimations
% Concentration values [Ideal model]
y_label = sprintf('Concentration (mol/m^3)');
title_name = sprintf('Concentration vs Time for H_2O - CO_2 mixture [Ideal State Model]');
legend_label = sprintf('Ideal State est');
% toSave = 'C_h2oco2_ideal.png';
text_label = sprintf('C(mol/m^3)[t] = %.4f + (%.4f)*exp(-%.4f*t)', x2IdealEst(1), x2IdealEst(2), x2IdealEst(3));
visualize1(C2Idealexp, C2Idealmodel, t2, title_name, y_label, legend_label, text_label, 20, 108, toSave);

% Concentration values [Ideal model]
y_label = sprintf('Concentration (mol/m^3)');
title_name = sprintf('Concentration vs Time for H_2O - CO_2 mixture [Real State Model]');
legend_label = sprintf('Real State est');
% toSave = 'C_h2oco2_real.png';
text_label = sprintf('C(mol/m^3)[t] = %.4f + (%.4f)*exp(-%.4f*t)', x2RealEst(1), x2RealEst(2), x2RealEst(3));
visualize1(C2Realexp, C2Realmodel, t2, title_name, y_label, legend_label, text_label, 20, 108, toSave);

% Pressure values
y_label = sprintf('P (bar)');
title_name = sprintf('Pressure vs Time for H_2O - CO_2 mixture [Both Models]');
% toSave = 'P_h2oco2.png';
visualize2(P2, (P2Idealest./bar_to_pascal), (P2Realest./bar_to_pascal), t2, title_name, y_label, toSave);


%% Model Adequacy checking (based on the P values) & displaying other results
% To display the total number of moles of CO2 being absorbed
n11 = x1IdealEst(2) * (V*mL_to_m3 / 2);
n12 = x1RealEst(2) * (V*mL_to_m3 / 2);
n21 = x2IdealEst(2) * (V*mL_to_m3 / 2);
n22 = x2RealEst(2) * (V*mL_to_m3 / 2);

% model 11 [part 1] ideal state model
RSq11 = rsqaured(P1*bar_to_pascal, P1Idealest);
fprintf('[Model 11]: EtOH-CO2 (Ideal Eq)\n');
fprintf('A = %.4f\nB = %.4f\nD = %.4f\n', x1IdealEst(1), x1IdealEst(2), x1IdealEst(3));
fprintf('R^2 = %.6f\n', RSq11);
fprintf('moles of CO2 being absorbed = %.8f\n\n', n11);

% model 12 [part 1] real state model
RSq12 = rsqaured(P1*bar_to_pascal, P1Realest);
fprintf('[Model 12]: EtOH-CO2 (Real Eq)\n');
fprintf('A = %.4f\nB = %.4f\nD = %.4f\n', x1RealEst(1), x1RealEst(2), x1RealEst(3));
fprintf('R^2 = %.6f\n', RSq12);
fprintf('moles of CO2 being absorbed = %.8f\n\n', n12);

% model 21 [part 2] ideal state model
RSq21 = rsqaured(P2*bar_to_pascal, P2Idealest);
fprintf('[Model 21]: H2O-CO2 (Ideal Eq)\n');
fprintf('A = %.4f\nB = %.4f\nD = %.4f\n', x2IdealEst(1), x2IdealEst(2), x2IdealEst(3));
fprintf('R^2 = %.6f\n', RSq21);
fprintf('moles of CO2 being absorbed = %.8f\n\n', n21);

% model 22 [part 2] real state model
RSq22 = rsqaured(P2*bar_to_pascal, P2Realest);
fprintf('[Model 22]: H2O-CO2 (Real Eq)\n');
fprintf('A = %.4f\nB = %.4f\nD = %.4f\n', x2RealEst(1), x2RealEst(2), x2RealEst(3));
fprintf('R^2 = %.6f\n', RSq22);
fprintf('moles of CO2 being absorbed = %.8f\n', n22);

fprintf('\n[Info] Units:\nA & B(mol/m^2), C (/s)\n');
%% Visualization Function 1 (plots for two models & the true one)
function visualize2(X, X_est1, X_est2, t, title_name, y_label, toSave)
    % X: the true values
    % X_est: the estimated values from different models
    % - X_est1: ideal state model
    % - X_est2: real state model
    % t: the x axis values

    figure();
    hold on;
    scatter(t, X, 50, 'blue', 'x', LineWidth=1.05);
    plot(t, X_est1, 'r-.', LineWidth=1.9);
    plot(t, X_est2, 'g-.', LineWidth=1.1);
    grid on;
    title(title_name);
    xlabel('Time (min)');
    ylabel(y_label);
    legend('Exp', 'Ideal State Est', 'Real State Est');

    if ~strcmp(toSave, 'none')
        saveas(gca, toSave);
    end
    hold off;
end

%% Visualization Function 2 (plots for one model & the true one)
function visualize1(X, X_est, t, title_name, y_label, legend_label, text_label, x1, y1, toSave)
    % X: the true values
    % X_est: the estimated values from a model
    % t: the x axis values

    figure();
    hold on;
    scatter(t, X, 50, 'blue', 'x', LineWidth=1.05);
    plot(t, X_est, 'r-.', LineWidth=1.1);
    grid on;
    title(title_name);
    xlabel('Time (min)');
    ylabel(y_label);
    legend('Exp', legend_label);

    if ~strcmp(text_label, 'none')
        text(x1, y1, text_label);
    end
    if ~strcmp(toSave, 'none')
        saveas(gca, toSave);
    end
    hold off;
end

%% Function for model adequacy checking
function R2 = rsqaured(X, Xest)
    % Find the R^2 of the model (model adequacy checking)
    % R^2 = 1 - (SS_res/SS_t) = SS_r / SS_t
    Xmean = mean(X);
    SS_t = sum((X - Xmean).^2);
    SS_r = sum((Xest - Xmean).^2);

    R2 = SS_r / SS_t;
end