%% CH3522 Experiment 4: Adsorption Kinetics
%  Author: Deepanjhan Das [CH22B020]
clear; close all;
format long;

%% Calibration Table

% wavelength (nm) recorded at
lambda = 663; 

% concentration (ppm)
conc = 0:1:5;

% absorbance 
absorbance = [0.0001, 0.1582, 0.3242, 0.5281, 0.7360, 0.9034];

% Performing linear regression with respect to a TLS objective assuming
% that there might error in both the parameters' calculations or observations
% conc (y) = m * absorbance (x) + c
[m, c, RSq] = TLSmodel(absorbance, conc);
conc_est = m * absorbance + c;


%% Visualize the data and the fit
figure();
hold on;
scatter(absorbance, conc, 40, 'blue', '+', LineWidth=1.1);
plot(absorbance, conc_est, 'r-', LineWidth=1);
xlabel('Absorbance (%)');
ylabel('Concentration (ppm)');
title('Calibration Table | Absorbance vs Concentration');
legend('Experimental Data', 'Fitted Data');
txt = sprintf('Model: (concentration) = %.6f * (absorbance) + %.6f\nR^2 = %.6f', m, c, RSq);
text(0.2, 4, txt);
grid on;
% save(gca, 'calibration.png');
hold off;


%% Function to fit the data to a TLS objective
function [m, c, RSq] = TLSmodel(x, y)
    covxy = cov(x,y,1);
    meanx = mean(x);
    meany = mean(y);

    m = ((covxy(2,2) - covxy(1,1)) + sqrt((covxy(2,2) - covxy(1,1))^2 + 4*(covxy(1,2)^2))) / (2*covxy(1,2));
    c = meany - (m * meanx);

    % Find the R^2 of the model (model adequacy checking)
    % R^2 = 1 - (SS_res/SS_t)
    y_est = m*x + c;
    SS_t = sum((y - meany).^2);
    SS_res = sum((y - y_est).^2);
    RSq = 1 - (SS_res / SS_t);
end