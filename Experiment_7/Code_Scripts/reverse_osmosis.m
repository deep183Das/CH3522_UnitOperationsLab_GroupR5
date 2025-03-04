%% CH3522 Experiment 7: Reverse Osmosis
%  Author: Deepanjhan Das [CH22B020]
clear; clc; close all;
format long;

%% Data 

% 1 ppm = 1 mg/L
% CF (ppm) : Feed Concentration
CF = [300, 300, 300, 400, 400, 400, 500, 500, 500]';
% CR (ppm) : Reject Concentration
CR = [508, 523, 554, 653, 680, 794, 805, 885, 1030]';
% CP (ppm) : Permeate Concentration
CP = [6, 5, 3, 5, 4, 6, 5, 6, 8]';

% PF : 1 (kgf/cm^2) = 0.967841 atm : External Pressure (atm)
PF = [40, 50, 60, 40, 50, 70, 45, 55, 75]';
% PR (atm)
PR = 10;
% PP (atm)
PP = 10;

% QR (L/s)
QR = [7.08, 6.85, 6.77, 6.89, 6.54, 4.73, 5.93, 5.36, 4.15]' ./ 1000;
% QP (L/s)
QP = [6.61, 5.66, 5.73, 6.93, 6.81, 7.22, 6.37, 7.21, 9.64]' ./ 1000;


%% Calculations

% rejection
Rej = 1 - (CP./CF);

% Non-ideality factor in π. (=2) for KCl as it dissociates into two ions.
phi = 2;

% Kelvin
T = 30 + 273.15;
% Universal Gas Constant (L.atm/mol.K)
R = 0.0821;
% Molecular weight of KCl (mg/mol)
MWkcl = 74.55 * (10^3);

% π : Osmotic pressure
Pi = @(C) ((phi .* C .* R .* T) ./ (MWkcl));

% π for feed
PiF = Pi(CF);
% π for reject
PiR = Pi(CR);
% π for permeate
PiP = Pi(CP);

P_FR = (PF + PR) ./ 2;
Pi_FR = (PiF + PiR) ./ 2;

del_P = (PF - PP) - (PiF - PiP);
delP_net = (P_FR - PP) - (Pi_FR - PiP);

% Jw : Flux of Water across membrane 
JwA = QP;  % multiplied with A. in L/s
% Js : Flux of Solute across membrane
JsA = CP .* QP .* ((10^3)./MWkcl); % multiplied with A. in mmol/s.

% kw : Mass Transfer coefficient for water
kwA = JwA ./ delP_net; % multiplied with A. L/s.atm
% ks : Mass Transfer coefficient for solute
delC = (CF - CP) ./ (MWkcl*10^-3);
ksA = JsA ./ delC; % multiplied with A. in mL/s.
ksA1e6 = ksA .* (10^6); % 10^-6 for each of them.

% r : Recovery
QF = QP + QR;
r = QP ./ QF;

% predicted CR 
CR_pred = CF .* ((1-((1-Rej).*r))./(1-r));

% relative % error
relError = abs((CR - CR_pred) ./ CR) .* 100;