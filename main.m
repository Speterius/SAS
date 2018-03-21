%% Init
clc; clear all; close all;
addpath('functions', 'ML_ae4304')

%% 0) GENERATE A/C + GUST + CONTROLLER STATE SPACE MODELS
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
% Graphically chosen gain for 0.5 damping ratio:
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

%% 1) STABILITY ANALYSIS
%%
% a) state space systems used: 
%display(uncont)
%display(cont)

% b) stab analysis:
figure; iopzmap(uncont(1:4,1)); title('Pitch damper OFF');
figure; iopzmap(cont(1:4,1)); title('Pitch damper ON');
% All longitudinal eigenmodes (short period + phugoid) are stable for both cases.

%% 2) TIME-DOMAIN ANALYSIS
%%
% Extend the state spaces with an 8th output: gamma_dot = theta_dot - alpha_dot
% // See extend_state_space.m //
uncont = extend_state_space(uncont);
cont = extend_state_space(cont);

%% Time Simulation Settings:
dt = 0.01;
T = 1000;

% set seed for consistent results:
seed = 32;

%Run time simulation with pitch-damper OFF:
% // See time_simulation.m //

[t, V_off, alpha_off, theta_off, qc_V_off, N_z_off] = time_simulation(uncont, dt, T, seed);
%Run time simulation with pitch-damper ON:
[~, V_on, alpha_on, theta_on, qc_V_on, N_z_on] = time_simulation(cont, dt, T, seed);

% Visualize:
figure; plot(t, V_off,      t, V_on);     legend('OFF', 'ON'); title('Velocity - V')
figure; plot(t, alpha_off,  t, alpha_on); legend('OFF', 'ON'); title('Angle of attack- \alpha')
figure; plot(t, theta_off,  t, theta_on); legend('OFF', 'ON'); title('Pitch - \theta')
figure; plot(t, qc_V_off,   t, qc_V_on);     legend('OFF', 'ON'); title('Pitch rate - q')
figure; plot(t, N_z_off,    t, N_z_on);   legend('OFF', 'ON'); title('Load factor - N_z')

%% 3) SPECTRAL ANALYSIS:
%% Analytical:


%% Experimental FFT:
%example51.m


%% Experimental PWELCH:



%% 4) VARIANCES:
%% Using analytical power spectra:


%% LYAPUNOV:
%ii) example72
