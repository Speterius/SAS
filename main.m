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
[cessna, ~, ~, ~, ~, V, c, sigmaug_V, sigmaag] = generate_state_space(sigma_wg, Lg);
%%
% Look at uncontrolled phugoid:
disp('Building aircraft model.')
uncont = cessna;
opt = stepDataOptions('StepAmplitude',-1);
figure; step(uncont(1:4, 1), opt); grid; title('Autopilot off')

% Add autopilot controller:
% Select K_theta from Root-locus with negative gains:
dK = 0.1;
%
%figure; rlocus(cessna(3, 1), 0:-dK:-50);title('Cessna pitch root locus')
% Graphically chosen gain for 0.5 damping ratio:
disp('Adding pitch damper.')
K_theta = -0.11703;
K = [0 0 K_theta 0 0 0 0];

% Construct controlled model:
cont = cessna;
cont.A = cessna.A-cessna.B(:,1)*K;
figure; step(cont(1:4, 1), opt); title('Autopilot on')

% Check phugoid damping:
[~, zeta] = damp(tf(cont(3, 1)));
zeta_phug = zeta(1);
display(zeta_phug)
%% 1) STABILITY ANALYSIS
%%
% a) state space systems used: 
% display(uncont(1:4, 1))
% display(cont(1:4, 1))

% b) stab analysis:
% aircraft is dynamically stable if both the phugoid and short period are stable ->
%       check poles are in LHP for beta to theta:

[~, ~, poles_uncont] = damp(uncont(3, 1));
[~, ~, poles_cont] = damp(cont(3, 1));

if any(real(poles_uncont)>0)
    disp('Uncontrolled system is unstable.')
else
    disp('Uncontrolled system is stable.')
end
if any(real(poles_cont)>0)
    disp('Controlled system is unstable.')
else
    disp('Controlled system is stable.')
end

% Show pole zero map for all input-output loops;
%figure; iopzmap(uncont(1:4,1)); title('Pitch damper OFF');
%figure; iopzmap(cont(1:4,1)); title('Pitch damper ON');

% All longitudinal eigenmodes (short period + phugoid) are stable for both cases.
%% 2) TIME-DOMAIN ANALYSIS
%%
% Extend the state spaces with an 8th output: N_z = V/g*(q - alpha_dot)
% // See extend_state_space.m //
uncont = extend_state_space(uncont, V, c);
cont = extend_state_space(cont, V, c);

% Time Simulation Settings:
dt = 0.01;
T = 1000;
seed = 32;

% Run time simulation with pitch-damper OFF:
% // See time_simulation.m //
disp('Running time simulation.')
[t, V_off, alpha_off, theta_off, qc_V_off, N_z_off] = time_simulation(uncont, dt, T, seed, V, c);
% Run time simulation with pitch-damper ON:
[~, V_on, alpha_on, theta_on, qc_V_on, N_z_on] = time_simulation(cont, dt, T, seed, V, c);

% Visualize Time plots for the 5 states:
figure; plot(t, V_off, t, V_on); title('Velocity - V'); ylabel('V [m/s]'); 
legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid;

figure; plot(t, alpha_off, t, alpha_on); title('Angle of attack - \alpha'); ylabel('\alpha [deg]'); 
legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid;

figure; plot(t, theta_off,  t, theta_on); title('Pitch angle - \theta'); ylabel('\theta [deg]'); 
legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid;
 
figure; plot(t, qc_V_off,   t, qc_V_on); title('Pitch rate - q'); ylabel('q [deg/s]'); 
legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid;

figure; plot(t, N_z_off,    t, N_z_on); title('Load factor - N_z'); ylabel('N_z [g units]');   
legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid; 

%% 3) SPECTRAL ANALYSIS:
% Subscripts:   - V/alpha/theta/q/N -for the 5 states
%               - off/on            -for pitch damper status
%               - a/e/p             -for analytical/experimental/pwelch
%
%                                   -> 30 PSDs
%% 3a) Analytical:
% // see analytic_psd.m //

% Frequency data:
w = logspace(-2,2,200);

% Uncontrolled aircraft:
S_V_off_a     = analytic_psd(uncont, 3, 1, w);
S_alpha_off_a = analytic_psd(uncont, 3, 2, w);
S_theta_off_a = analytic_psd(uncont, 3, 3, w);
S_q_off_a     = analytic_psd(uncont, 3, 4, w);
S_N_off_a     = analytic_psd(uncont, 3, 8, w);

% Controlled aircraft
S_V_on_a     = analytic_psd(cont, 3, 1, w);
S_alpha_on_a = analytic_psd(cont, 3, 2, w);
S_theta_on_a = analytic_psd(cont, 3, 3, w);
S_q_on_a     = analytic_psd(cont, 3, 4, w);
S_N_on_a     = analytic_psd(cont, 3, 8, w);

%% 3b) 3c) Experimental using FFT and pwelch:

dt = 0.1;
T = 10000;
NFFT = 2048;
disp('Running PSD calculation.')
[omega, S_V_off_e, S_alpha_off_e, S_theta_off_e, S_q_off_e, S_N_off_e, PW_uncont, PW_w_uncont] = experi_psd(uncont, dt, T+dt, seed, NFFT);
[~,     S_V_on_e,  S_alpha_on_e,  S_theta_on_e,  S_q_on_e,  S_N_on_e, PW_cont, PW_w_cont]  = experi_psd(cont, dt, T+dt, seed, NFFT);

% Collect data from PWELCH output matrix:
S_V_off_p       = PW_uncont(:, 1);
S_alpha_off_p   = PW_uncont(:, 2);
S_theta_off_p   = PW_uncont(:, 3);
S_q_off_p       = PW_uncont(:, 4);
S_N_off_p       = PW_uncont(:, 5);

S_V_on_p        = PW_cont(:, 1);
S_alpha_on_p    = PW_cont(:, 2);
S_theta_on_p    = PW_cont(:, 3);
S_q_on_p        = PW_cont(:, 4);
S_N_on_p        = PW_cont(:, 5);

% rename frequency axes:
w_a = w;
w_e = omega;
w_p = 2*pi*PW_w_cont;

%% 3d) Visualize PSDS plots:
plot_psds
plot_psds_subplots

%% 4) VARIANCES:
%% 4a) Using analytical power spectra:


%% 4b) Lyapunov:
%ii) example72
