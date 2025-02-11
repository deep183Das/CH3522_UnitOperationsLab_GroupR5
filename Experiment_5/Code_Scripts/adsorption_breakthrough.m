%% CH3522 Experiment 5: Adsorption Breakthrough
%  Author: Deepanjhan Das [CH22B020]
clear; close all;
format long;

%% Data of RH & Time

% time in minute
t = 0:1:54;
% relative humidity (%)
rh = [74.1, 59.3, 56.7, 58.4, 60.5, 62.2, 64.5, 66.1, 66.1, 66.8, 67.4, ...
    67.9, 68.8, 69.6, 69.7, 70.1, 70.3, 70.7, 70.9, 71.4, 71.8, 71.9, ...
    72.1, 71.9, 72.1, 72.1, 72.3, 72.4, 72.6, 72.8, 72.8, 72.9, 72.8, ...
    72.9, 73.0, 73.2, 73.2, 73.3, 73.3, 73.2, 73.2, 73.0, 72.9, 72.9, ...
    72.9, 72.8, 72.7, 72.7, 72.6, 72.6, 72.4, 72.2, 71.7, 71.4, 69.5];

% A. The (RH vs Time) plot
figure();
hold on;
scatter(t, rh, 40, 'blue', '+', LineWidth=1.1);
grid on;
xlabel('Time (minutes)');
ylabel('Relative Humidity (%)');
title('Trend on Relative Humidity with Time');
legend('Experimentally Observed Data');
hold off;


%% Approximation of C_out / C_initial from the above RH data.
%  Each of the C terms is denoting mass concentrations.
ratio_of_Co_to_Ci = rh ./ (rh(1));
VF = 3.37078 * 10^-6;  % in m^3/s

% B. The (Co/Ci vs VF*t) plot where VF is the volumetric flow rate of air
figure();
hold on;
scatter(t.*60, ratio_of_Co_to_Ci, 40, '+', 'blue', LineWidth=1.1);
xlabel('Time (minutes)');
ylabel('C_t/C_0');
title('Relative mass concentrations trend with Time');
legend('Experimentally Observed Data');
grid on;
hold off;


%% Defining paramters and constants for ROSEN Model. 
%  To fit the model, we need to perform a non-linear optimization using an
%  initial value of the overall MTC whose fitted estimate will be obtained.

T = 30 + 273.15;   % experiment temperature (K)
Dai = 2.3 * 10^-5; % diffusion coefficient of water vapour (m^2/s)
R = 0.0051;        % column radius (1/2 of the diameter in m)
Uz = VF/(pi*R^2);  % superficial velocity of the vapour (m/s)
Z = 0.057;         % packing height (m)
m0 = 100;          % numerical constant in Rosen Model (initial value)
                   % dp : silica gel particle size (63-200 nm, 70-230 mesh)
rho_s = 2200;      % density of silica gel (Kg/m^3)
Ka3i = 20;         % linear adsorption constant (m^3/Kg)

% bed length parameter (without m)
X = (3*Dai*Ka3i*Z*rho_s)/(10*Uz*R^2);
Y = @(t) ((2.*Dai)./(R.^2)).*(t.*60 - (Z./Uz));
Yt = Y(t);
% nu/K is the film resistance parameter
nu = (Dai*Ka3i*rho_s*10)/R;

% initial estimate of K (m/s)
K0 = 0.01;  
x0(1) = K0;
x0(2) = m0;
x10(1) = 1;
x10(2) = 1;
x10(3) = 1;
x10(4) = 1;
options = optimoptions(@lsqnonlin, 'Algorithm','levenberg-marquardt', 'Display','iter');
x = lsqnonlin(@(x)ObjectiveFunc(x,X,Yt(2:end),nu,ratio_of_Co_to_Ci(2:end)), x0, 1e-4,[], options);
x1 = lsqnonlin(@(x1)ObjectiveFunc1(x1,t(2:end)*60,ratio_of_Co_to_Ci(2:end)), x10, [],[], options);

% [Result 1] Estimation of Overall MTC using breakthrough curve
fprintf('Overall MTC (K_Ai) is = %.6f\n', x(1));

% % find the fitted values
% est_ratio_of_Co_to_Ci = 0.5.*(1 + erf((((3.*Y(t))./(2.*X))-1)./(2.*sqrt((nu)./(X.*x(1))))));
% % find the R2 of the model
% R2 = RSq(ratio_of_Co_to_Ci, est_ratio_of_Co_to_Ci);
% fprintf('R^2 of the model = %.6f\n', R2);

t_new = linspace(0,52,201)*60; % seconds
% new_ratio_of_Co_to_Ci = 0.5.*(1 + erf((((3.*Y(t_new))./(2.*X))-1)./(2.*sqrt((nu)./(X.*x(1))))));
% new_ratio = x1(1) ./ (x1(2) + exp(-t_new.*x1(3)));
new_ratio = x1(1) + (x1(2).*erf((t_new-x1(3))/(x1(4))));

% The adsorption breakthrough curve (use the fitted data too) 
figure();
hold on;
scatter(t*60, ratio_of_Co_to_Ci, 40, 'blue', '+', LineWidth=1.1);
plot(t_new, new_ratio, 'r-', LineWidth=1.2);
grid on;
legend('Exp Data', 'Fitted Data');
xlabel('Time (second)');
ylabel('C_t/C_0');
title('C_t/C_0 vs time curve');
legend('Experimental Data');
hold off;

% [Result 2] Estimation of maximum (equilibrium) adsorption 
% capacity of the silica gel bed with water.

RH = rh(2:end) ./ 100;
P_sat = 4239.651; % Pa
Rv = 461.5;  % J/Kg.K
rho_v = (RH .* P_sat) ./ (Rv * T);

y_integral = VF .* rho_v .* (1 - ratio_of_Co_to_Ci(2:end));
x_integral = t(2:end) .* 60;
% use trapezoid rule
mass = trapz(x_integral, y_integral);
fprintf('Mass Adsorbed = %.16f Kg\n', mass);

% %% Function to compute the R^2 of the model
% function R2 = RSq(y, y_est)
%     % Find the R^2 of the model (model adequacy checking)
%     meany = mean(y);
% 
%     % R^2 = 1 - (SS_res/SS_t)
%     SS_t = sum((y - meany).^2);
%     SS_res = sum((y - y_est).^2);
%     R2 = 1 - (SS_res / SS_t);
% end