%% Init
clc; clear all; close all;
addpath('functions', 'ML_ae4304')

T = 100;
dt = 0.1;

%% 0) Generate aircraft + turbulence state space models
%%
% Turbulence scale length and intensity::
Lg          = 150;      % [m]
sigma_wg    = 1;        % [m/s]

% Generate A and B matrices using aircraft c&s derivatives provided:
% (see generate_state_space.() )
[cessna, ~, ~, ~, ~, V, c, sigmaug_V, sigmaag] = generate_state_space(Lg, sigma_wg);


%%
% Look at uncontrolled phugoid:
uncont = cessna(:, 1);
opt = stepDataOptions('StepAmplitude',-1);
figure; step(uncont(1:4, 1), opt); grid; title('Autopilot off')

% Add autopilot controller:
cont = cessna;
K_theta = -0.3;
K = [0 0 K_theta 0 0 0 0];
cessna_cont.A = cessna.A-cessna.B(:,1)*K;

cont   = cont(1:4, 1);
figure; step(cont, opt); title('Autopilot on')

%% TUNE GAIN BY POLE PLACEMENT






%% 1) Stability analysis
%%



%% 2) Time-Domain Analysis
%%



%% 3) Spectral Analysis
%%