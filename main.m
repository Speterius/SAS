%% Init
clc; clear all; close all;
addpath('functions', 'ML_ae4304', 'plotting_scripts')
str = input('Show all 30 figures?[y/n]: ','s');
if isempty(str)
    str = 'Y';
end
if str == 'y' || str == 'Y'
    visualize = true;
elseif str == 'n' || str == 'N'
    visualize = false;
end

%% 0) GENERATE A/C + GUST + CONTROLLER STATE SPACE MODELS
%%
% Turbulence scale length and intensity::
Lg          = 150;      % [m]
sigma_wg    = 1;        % [m/s]

% Generate A and B matrices using aircraft c&s derivatives provided:
% // see generate_state_space.m //
[cessna, ~, ~, ~, ~, V, c, sigmaug_V, sigmaag] = generate_state_space(sigma_wg, Lg);

% Uncontrolled phugoid step response:
disp('Building aircraft model.')
uncont = cessna;
opt = stepDataOptions('StepAmplitude',-1);

if visualize
    figure; step(uncont(3, 1), opt); grid; title('Autopilot off')
end
    
% Add autopilot controller:
% Select K_theta from Root-locus with negative gains:
dK = 0.1;
if visualize
    figure; rlocus(cessna(3, 1), 0:-dK:-50);title('Cessna pitch root locus')
end

% Graphically chosen gain for 0.5 damping ratio:
disp('Adding pitch damper.')
K_theta = -0.11703;
K = [0 0 K_theta 0 0 0 0];

% Construct controlled model and check step response:
cont = cessna;
cont.A = cessna.A-cessna.B(:,1)*K;

if visualize
    figure; step(cont(3, 1), opt); title('Autopilot on'); grid; 
end

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

[~, ~, poles_uncont] = damp(uncont(:, :));
[~, ~, poles_cont] = damp(cont(:, :));

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

%% Time Simulation Settings:
dt = 0.1;
T = 5000;
seed = 32;

% Run time simulation with pitch-damper OFF:
% // See time_simulation.m //
disp('Running time simulation.')
[t, V_off, alpha_off, theta_off, qc_V_off, N_z_off] = time_simulation(uncont, dt, T, seed);
% Run time simulation with pitch-damper ON:
[~, V_on, alpha_on, theta_on, qc_V_on, N_z_on] = time_simulation(cont, dt, T, seed);

% Visualize Time plots for the 5 states:
if visualize
    figure; plot(t, V_off, t, V_on); title('Velocity - V'); ylabel('u/V [-]'); 
    legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid;

    figure; plot(t, alpha_off, t, alpha_on); title('Angle of attack - \alpha'); ylabel('\alpha [rad]'); 
    legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid;

    figure; plot(t, theta_off,  t, theta_on); title('Pitch angle - \theta'); ylabel('\theta [rad]'); 
    legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid;

    figure; plot(t, qc_V_off,   t, qc_V_on); title('Pitch rate - q'); ylabel('qc/V [-]'); 
    legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid;

    figure; plot(t, N_z_off,    t, N_z_on); title('Load factor - N_z'); ylabel('N_z [g units]');   
    legend('Pitch damper OFF', 'Pitch damper ON'); xlabel('t [s]'); grid;
end

%% 3) SPECTRAL ANALYSIS:
% Subscripts:   - V/alpha/theta/q/N -for the 5 states
%               - off/on            -for pitch damper status
%               - a/e/p             -for analytical/experimental/pwelch
%
%                                   -> 30 PSDs, lots of figures
%% 3a) Analytical:
% // see analytic_psd.m //
% analytic_psd() - input_nr = 3 for vertical gust:

% Frequency data:
Nomega = 2000;
w_a = logspace(-2,2,Nomega);

% Uncontrolled aircraft:
S_V_off_a     = analytic_psd(uncont, 3, 1, w_a);
S_alpha_off_a = analytic_psd(uncont, 3, 2, w_a);
S_theta_off_a = analytic_psd(uncont, 3, 3, w_a);
S_q_off_a     = analytic_psd(uncont, 3, 4, w_a);
S_N_off_a     = analytic_psd(uncont, 3, 8, w_a);

% Controlled aircraft
S_V_on_a     = analytic_psd(cont, 3, 1, w_a);
S_alpha_on_a = analytic_psd(cont, 3, 2, w_a);
S_theta_on_a = analytic_psd(cont, 3, 3, w_a);
S_q_on_a     = analytic_psd(cont, 3, 4, w_a);
S_N_on_a     = analytic_psd(cont, 3, 8, w_a);

S_xx_off = [S_V_off_a, S_alpha_off_a, S_theta_off_a, S_q_off_a, S_N_off_a];
S_xx_on = [S_V_on_a, S_alpha_on_a, S_theta_on_a, S_q_on_a, S_N_on_a];

%% 3b) Experimental using FFT:

% Settings also affect pwelch:
seed = 32;
dt = 0.1;
T = 1000;
disp('Running DFT calculations.')
[omega, S_V_off_e, S_alpha_off_e, S_theta_off_e, S_q_off_e, S_N_off_e] = experi_psd(uncont,dt, T+dt, seed);
[~,     S_V_on_e,  S_alpha_on_e,  S_theta_on_e,  S_q_on_e,  S_N_on_e]  = experi_psd(cont,  dt, T+dt, seed);
w_e = omega;

%% 3c) PWELCH Routine:
% // based on example 4.5 // 
disp('Running PWELCH.')

% Settings:
N = T/dt;
window = 500;
noverlap = 25;
fs = 1/dt;

% Get time traces from lsim:
[t_array, v_off, alpha_off, theta_off, q_off, N_off] = time_simulation(uncont, dt, T, seed);
[~, v_on, alpha_on, theta_on, qc_V_on, N_z_on] = time_simulation(cont, dt, T, seed);

% Execute pwelch.m
[S_V_off_p,    f_p] = pwelch(v_off,     window, noverlap, N, fs, 'onesided'); S_V_off_p=S_V_off_p/2; S_V_off_p = S_V_off_p(2:N/2+1);
[S_alpha_off_p, ~ ] = pwelch(alpha_off, window, noverlap, N, fs, 'onesided'); S_alpha_off_p=S_alpha_off_p/2; S_alpha_off_p = S_alpha_off_p(2:N/2+1);
[S_theta_off_p, ~ ] = pwelch(theta_off, window, noverlap, N, fs, 'onesided'); S_theta_off_p=S_theta_off_p/2; S_theta_off_p = S_theta_off_p(2:N/2+1);
[S_q_off_p,     ~ ] = pwelch(q_off,     window, noverlap, N, fs, 'onesided'); S_q_off_p=S_q_off_p/2; S_q_off_p = S_q_off_p(2:N/2+1);
[S_N_off_p,     ~ ] = pwelch(N_off,     window, noverlap, N, fs, 'onesided'); S_N_off_p=S_N_off_p/2; S_N_off_p = S_N_off_p(2:N/2+1);

[S_V_on_p,      ~ ] = pwelch(v_on,      window, noverlap, N, fs, 'onesided'); S_V_on_p=S_V_on_p/2; S_V_on_p = S_V_on_p(2:N/2+1);
[S_alpha_on_p,  ~ ] = pwelch(alpha_on,  window, noverlap, N, fs, 'onesided'); S_alpha_on_p=S_alpha_on_p/2; S_alpha_on_p = S_alpha_on_p(2:N/2+1);
[S_theta_on_p,  ~ ] = pwelch(theta_on,  window, noverlap, N, fs, 'onesided'); S_theta_on_p=S_theta_on_p/2; S_theta_on_p = S_theta_on_p(2:N/2+1);
[S_q_on_p,      ~ ] = pwelch(qc_V_on,   window, noverlap, N, fs, 'onesided'); S_q_on_p=S_q_on_p/2; S_q_on_p = S_q_on_p(2:N/2+1);
[S_N_on_p,      ~ ] = pwelch(N_z_on,    window, noverlap, N, fs, 'onesided'); S_N_on_p=S_N_on_p/2; S_N_on_p = S_N_on_p(2:N/2+1);

% Handle frequency axis:
f_p = f_p(2:N/2+1);
w_p = 2*pi*f_p;

% SHOW ALL separate 15 figures for the 30 PSDS:
disp('done')
if visualize
    plot_psd_analytic
    plot_psd_fft
    plot_psd_welch
end

%% 3d) Show figures for COPMARISON of methods:
if visualize
    plot_psds_subplots
end

%% 4) VARIANCES:
disp('Calculating variances.')
%% 4a) Using analytical power spectra:
% // based on example74a.m //
% Sxx = [S_v S_alpha S_theta S_q S_N] 

n_states = 5;
var_off_a = zeros(1,n_states);
var_on_a = zeros(1,n_states);

% Integrate power spectral density for all 5 states:
for i=1:Nomega-1
    for j=1:n_states
        var_off_a(j)=var_off_a(j)+(w_a(i+1)-w_a(i))*S_xx_off(i,j);
        var_on_a(j)=var_on_a(j)+(w_a(i+1)-w_a(i))*S_xx_on(i,j);
    end
end

var_off_a=var_off_a/pi;
var_on_a = var_on_a/pi;

%% 4b) Lyapunov:
% // based on example72.m //

% Intensity is 1:
Wc = 1.0;
% Pitch damper OFF: uncont
A = uncont.A;
B = uncont.B(:,3);
C = uncont.C;
D = uncont.D;

L_off = lyap(A,B*Wc*B');


% LOAD FACTOR VARIANCE:
C_N = C(8,:);
D_N = D(8,:);

% Only interested in 3rd input channel:
W_N = [0, 0, 0 ; 
     0, 0, 0; 
     0, 0, 1];

var_NZ_off = C_N*L_off*C_N' + D_N*W_N*D_N';

% Select diagonals (1:4) of L + the load factor diagonal (8) of P:
lyap_var_off = [diag(L_off(1:4, 1:4))', var_NZ_off];

% Pitch damper ON: cont
A = cont.A;
B = cont.B(:,3);
C = cont.C;
D = cont.D;

L_on = lyap(A,B*Wc*B');

% LOAD FACTOR VARIANCE:
C_N = C(8,:);
D_N = D(8,:);

% Only interested in 3rd input channel:
W_N = [0, 0, 0 ; 
     0, 0, 0; 
     0, 0, 1];

var_NZ_on = C_N*L_on*C_N' + D_N*W_N*D_N';

% Select diagonals (1:4) of L + the load factor diagonal (8) of P:
lyap_var_on = [diag(L_on(1:4, 1:4))', var_NZ_on];

%% 4c) USING var.m:
% Get time traces:
dt = 0.01;
T = 100000;
seed = 32;
% Run time simulation with pitch-damper OFF:
[t, V_off, alpha_off, theta_off, qc_V_off, N_z_off] = time_simulation(uncont, dt, T, seed);
[~, V_on, alpha_on, theta_on, qc_V_on, N_z_on] = time_simulation(cont, dt, T, seed);

TRACE_off    = [V_off, alpha_off, theta_off, qc_V_off, N_z_off];
TRACE_on  = [V_on, alpha_on, theta_on, qc_V_on, N_z_on];

var_varm_off  = var(TRACE_off, 1);
var_varm_on   = var(TRACE_on, 1);

% Create table with variances:

% Pitch damper off:
VAR = [var_off_a; lyap_var_off; var_varm_off];
SIGMA_off = table(VAR(:, 1), VAR(:, 2), VAR(:, 3), VAR(:, 4), VAR(:, 5),...
    'VariableNames', {'sigma2_V','sigma2_alpha','sigma2_theta','sigma2_q','sigma2_Nz'},...
    'RowNames', {'From Analyitical PSD','From LYAP','From var.m'});

% Pitch damper on:
VAR = [var_on_a; lyap_var_on; var_varm_on];
SIGMA_on = table(VAR(:, 1), VAR(:, 2), VAR(:, 3), VAR(:, 4), VAR(:, 5),...
    'VariableNames', {'sigma2_V','sigma2_alpha','sigma2_theta','sigma2_q','sigma2_Nz'},...
    'RowNames', {'From Analyitical PSD','From LYAP','From var.m'});

display(SIGMA_off)
display(SIGMA_on)

%% Run for high N_it to get true variance estimate:
% // see variance_ensemble.m //

N_it    = 100;          % [-]
T_real  = 500;          % [s]
dt      = 0.1;          % [s]
print_percentage_done = false;

var_ensemble_on  = variance_ensemble(uncont, N_it, T_real, dt, print_percentage_done); display(var_ensemble_on);
var_ensemble_off = variance_ensemble(cont, N_it, T_real, dt, print_percentage_done); display(var_ensemble_off);

%% End session
if visualize
    disp('  Press any button to close all windows.')
    pause
end
close all;

%% Run plotting scripts separately to look at PSD figures:
%%
% plot_psd_analytic
%%
% plot_psd_fft
%%
% plot_psd_welch
%%
% plot_psds_subplots




