%% Init
clc; clear all; close all;
addpath('functions', 'ML_ae4304')

%% 0) Generate aircraft + turbulence state space models
%%
% Turbulence scale length and intensity::
Lg          = 150;      % [m]
sigma_wg    = 1;        % [m/s]

% Generate A and B matrices using aircraft c&s derivatives provided:
% // see generate_state_space.m //
[cessna, ~, ~, ~, ~, V, c, sigmaug_V, sigmaag] = generate_state_space(Lg, sigma_wg);

% Look at uncontrolled phugoid:
uncont = cessna;
opt = stepDataOptions('StepAmplitude',-1);
figure; step(uncont(1:4, 1), opt); grid; title('Autopilot off')

% Add autopilot controller:
% Select K_theta from Root-locus with negative gains:
dK = 0.1;
figure; rlocus(cessna(3, 1), 0:-dK:-50);title('Cessna pitch root locus')
K_theta = -0.4617;
K = [0 0 K_theta 0 0 0 0];

% Construct controlled model:
cont = cessna;
cont.A = cessna.A-cessna.B(:,1)*K;
figure; step(cont(1:4, 1), opt); title('Autopilot on')

% Check phugoid damping:
[~, zeta] = damp(cont(1:3, 1));
zeta_phug = zeta(3);
display(zeta_phug)

%% 1) Stability analysis
%%
% a) state space systems used: 
%display(uncont)
%display(cont)

% b) stab analysis:
figure; iopzmap(uncont(1:4,1)); title('Pitch damper OFF');
figure; iopzmap(cont(1:4,1)); title('Pitch damper ON');
% All longitudinal eigenmodes (short period + phugoid) are stable for both cases.

%% 2) Time-Domain Analysis
%%
% Extend the state spaces with an 8th output: gamma_dot = theta_dot - alpha_dot
% // See extend_state_space.m //
uncont = extend_state_space(uncont);
cont = extend_state_space(cont);

%% TIME SIMULATION SETTINGS:
dt = 0.01;
T = 100;

% set seed to get same results:
seed = 32;

%Pitch damper OFF:
[t, V_off, alpha_off, theta_off, q_off, N_z_off] = time_simulation(uncont, dt, T, seed);

% Pitch damper ON:
[~, V, alpha, theta, qc_V, N_z] = time_simulation(cont, dt, T, seed);

% VISUALIZE:
figure; plot(t, N_z_off, t, N_z); legend('OFF', 'ON')

%% 3) Spectral Analysis
%%

%example51.m



