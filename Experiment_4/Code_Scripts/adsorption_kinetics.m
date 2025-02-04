%% CH3522 Experiment 4: Adsorption Kinetics
%  Author: Deepanjhan Das [CH22B020]
clear; close all;
format long;

%% The model obtained from the calibration is as follows:
%  conc (ppm) = (5.422738524421927) * absorbance + 0.104957151713649
conc_est = @(absorbance) (5.422738524421927 .* absorbance + 0.104957151713649);

% observed time data (seconds)
time = [0, 5:5:25, 35:10:75, 90] * 60;

% obtained absorbance values (%)
absorbance_val = [0.1727, 0.1331, 0.5568, 0.5492, 0.4509, 0.2783, ...
    0.2439, 0.2680, 0.1412, 0.1629, 0.1200];

% dilution level to obtain the corrected concentration (ppm)
dil_level = [ones(1,2)*100, ones(1,9)*20];


%% The model is used to obtain the concentration of each sample:
concentration_est = [100.0, (conc_est(absorbance_val) .* dil_level)];
% Now using the model obtained from the calibration data, we get the
% estimated concentration (ppm) which correspond to the diluted version of
% the actual samples of concern and therefore we need to multiply them with
% the specific dilution level factors to obtain the actual or corrected
% estimations of the samples' concentrations at those specific times. These
% will give us the concentration of the dye in the system that was getting
% adsorbed in the activated charcoal used.


%% Non-linear optimization (OLS objective)
x0(1) = 100;
x0(2) = 0.02;
x0(3) = 60;
y = concentration_est;
t = time;

options = optimoptions(@lsqnonlin, 'Algorithm','levenberg-marquardt', 'Display','iter');
x = lsqnonlin(@(x)ConcFunc(x,y,t), x0, [],[], options);

% extract the important parameters from the fitted model
C0 = 100; % ppm
m = 4000; % mg activated charcoal (adsorbent)
V = 1.0;  % liter

k = (((V*C0) / x(3)) - V) / m;  % L/mg
Kca = x(2) / (1 + (V/(m*k)));    % /s
fprintf('\n[Info]\nEquilibrium constant: k = %.6f L/mg (or /ppm)\n', k);
fprintf('Liquid mass transfer coefficient: (Kc*a) = %.10f /s\n', Kca);

% estimate using the fitted model
y_est = x(1).*exp(-x(2).*t) + x(3);
R2 = RSq(y, y_est);
fprintf('Model R^2 = %.6f\n', R2);


%% Visualize the trend in the concentration change with time
figure();
hold on;
scatter(time, concentration_est, 40, 'blue', '+', LineWidth=1.1);
% plot the fitted model
t_new = linspace(t(1), t(end), 1001);
y_new = x(1).*exp(-x(2).*t_new) + x(3);
plot(t_new, y_new, 'r-', LineWidth=1.2);
legend('Experimental Val', 'Fitted Data');
title('Trend of change in concentration of the dye in the system');
xlabel('Time (seconds)');
ylabel('Concentration (ppm)');
txt = sprintf('Model: (concentration) = (%.6f)*exp(-%.6f * time) + %.6f\nModel Adequacy: R^2 = %.6f', x(1), x(2), x(3), R2);
text(1000, 80, txt);
grid on;
hold off;


%% Function to compute the R^2 of the model
function R2 = RSq(y, y_est)
    % Find the R^2 of the model (model adequacy checking)
    meany = mean(y);

    % R^2 = 1 - (SS_res/SS_t)
    SS_t = sum((y - meany).^2);
    SS_res = sum((y - y_est).^2);
    R2 = 1 - (SS_res / SS_t);
end