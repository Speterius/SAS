%% Init
clc; clear all; close all;
addpath('functions', 'ML_ae4304', 'plotting_scripts')

str = input('Show all figures? [y/n]: ','s');
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
%
% Look at uncontrolled phugoid:
disp('Building aircraft model.')
uncont = cessna;
opt = stepDataOptions('StepAmplitude',-1);

if visualize
    figure; step(uncont(3, 1), opt); grid; title('Autopilot off')
end
    
% Add autopilot controller:
% Select K_theta from Root-locus with negative gains:
dK = 0.1;
%
if visualize
    figure; rlocus(cessna(3, 1), 0:-dK:-50);title('Cessna pitch root locus')
end

% Graphically chosen gain for 0.5 damping ratio:
disp('Adding pitch damper.')
K_theta = -0.11703;
K = [0 0 K_theta 0 0 0 0];

% Construct controlled model:
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
[t, V_off, alpha_off, theta_off, qc_V_off, N_z_off] = time_simulation(uncont, dt, T, seed, V, c);
% Run time simulation with pitch-damper ON:
[~, V_on, alpha_on, theta_on, qc_V_on, N_z_on] = time_simulation(cont, dt, T, seed, V, c);

% Visualize Time plots for the 5 states:
if visualize
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
end

%% 3) SPECTRAL ANALYSIS:
% Subscripts:   - V/alpha/theta/q/N -for the 5 states
%               - off/on            -for pitch damper status
%               - a/e/p             -for analytical/experimental/pwelch
%
%                                   -> 30 PSDs
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

% ANALYTIC PLOTS:
if visualize
    plot_psd_analytic
end

%% 3b) Experimental using FFT:
seed = 32;
dt = 0.01;
T = 50000;
disp('Running DFT calculations.')
%[omega, S_V_off_e, S_alpha_off_e, S_theta_off_e, S_q_off_e, S_N_off_e] = experi_psd(uncont,dt, T+dt, seed);
%[~,     S_V_on_e,  S_alpha_on_e,  S_theta_on_e,  S_q_on_e,  S_N_on_e]  = experi_psd(cont,  dt, T+dt, seed);


% 3c) PWELCH Routine:
disp('Running PWELCH.')
% settings:
N = T/dt;
window = 0.1*N;
noverlap = 50;
NFFT = 2048;
fs = 1/dt;

fres = fs/N;

% Get time traces from lsim:
[t_array, v_off, alpha_off, theta_off, q_off, N_off] = time_simulation(uncont, dt, T, seed, V, c);
[~, v_on, alpha_on, theta_on, qc_V_on, N_z_on] = time_simulation(cont, dt, T, seed, V, c);

% Execute pwelch.m
%%
[S_V_off_p,    f_p] = pwelch(v_off,     window, noverlap, N, fs, 'onesided'); S_V_off_p=S_V_off_p/2; S_V_off_p = S_V_off_p(2:N/2+1);
[S_alpha_off_p, ~ ] = pwelch(alpha_off, window, noverlap, N, fs, 'onesided'); S_alpha_off_p=S_alpha_off_p/2; S_alpha_off_p = S_alpha_off_p(2:N/2+1);
[S_theta_off_p, ~ ] = pwelch(theta_off, window, noverlap, N, fs, 'onesided'); 
[S_q_off_p,     ~ ] = pwelch(q_off,     window, noverlap, N, fs, 'onesided'); 
[S_N_off_p,     ~ ] = pwelch(N_off,     window, noverlap, N, fs, 'onesided'); 

[S_V_on_p,      ~ ] = pwelch(v_on,      window, noverlap, N, fs, 'onesided'); S_V_on_p=S_V_on_p/2; S_V_on_p = S_V_on_p(2:N/2+1);
[S_alpha_on_p,  ~ ] = pwelch(alpha_on,  window, noverlap, N, fs, 'onesided');
[S_theta_on_p,  ~ ] = pwelch(theta_on,  window, noverlap, N, fs, 'onesided');
[S_q_on_p,      ~ ] = pwelch(qc_V_on,   window, noverlap, N, fs, 'onesided');
[S_N_on_p,      ~ ] = pwelch(N_z_on,    window, noverlap, N, fs, 'onesided');

% rename frequency axes:
%w_e = omega;
f_p = f_p(2:N/2+1);
w_p = 2*pi*f_p;
%
if visualize
    %plot_psd_fft
    %plot_psd_welch
end

%%
% pwelch Velocity
figure('rend','painters','pos',[10 10 900 600]); 
loglog(w_p, S_V_off_p, w_p, S_V_on_p, 'Linewidth', 1);
title('Welch PSD Velocity PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{VV}'); xlabel('log w'); grid;

%% 3d) Visualize PSDS plots:
if visualize
    plot_psds
    plot_psds_subplots
end

%% 4) VARIANCES:
% Clear memory except the state space models and the analyitical power spectra:
clearvars -except cont uncont w_a S_xx_off S_xx_on Nomega T dt V c visualize
%% 4a) Using analytical power spectra:
disp('Calculating variances.')
% // based on example74a.m //
% Sxx = [S_v S_alpha S_theta S_q S_N] 

n_states = 5;

var_off = zeros(1,n_states);
var_on = zeros(1,n_states);

for i=1:Nomega-1
    for j=1:n_states
        var_off(j)=var_off(j)+(w_a(i+1)-w_a(i))*S_xx_off(i,j);
        var_on(j)=var_on(j)+(w_a(i+1)-w_a(i))*S_xx_on(i,j);
    end
end

var_off=var_off/pi;
var_on = var_on/pi;

% 4b) Lyapunov:
% // based on example72.m //
Wc = 1.0;
% Pitch damper OFF: uncont
A = uncont.A;
B = uncont.B(:,3);
C = uncont.C;
D = uncont.D;

L_off = lyap(A,B*Wc*B');
% Output covariance matrix:
P_off = C*L_off*C' + D*Wc*D';

% Select diagonals (1:4) of L + the load factor diagonal (8) of P:
lyap_var_off = [diag(L_off(1:4, 1:4))', P_off(8, 8)];

% Pitch damper ON: cont
A = cont.A;
B = cont.B(:,3);
C = cont.C;
D = cont.D;

L_on = lyap(A,B*Wc*B');
% Output covariance matrix:
P_on = C*L_on*C' + D*Wc*D';

% Select diagonals (1:4) of L + the load factor diagonal (8) of P:
lyap_var_on = [diag(L_on(1:4, 1:4))', P_on(8, 8)];

%% 4c) USING var.m:
% Get time traces:
dt = 0.01;
T = 100000;
seed = 32;
% Run time simulation with pitch-damper OFF:
[t, V_off, alpha_off, theta_off, qc_V_off, N_z_off] = time_simulation(uncont, dt, T, seed, V, c);
[~, V_on, alpha_on, theta_on, qc_V_on, N_z_on] = time_simulation(cont, dt, T, seed, V, c);

TRACE_off    = [V_off, alpha_off, theta_off, qc_V_off, N_z_off];
TRACE_on  = [V_on, alpha_on, theta_on, qc_V_on, N_z_on];

weight = 1;
V_off  = var(TRACE_off, weight);
V_on   = var(TRACE_on, weight);


% Create table with variances:

% Pitch damper off:
VAR = [var_off; lyap_var_off; V_off];
SIGMA_off = table(VAR(:, 1), VAR(:, 2), VAR(:, 3), VAR(:, 4), VAR(:, 5),...
    'VariableNames', {'sigma_V','sigma_alpha','sigma_theta','sigma_q','sigma_Nz'},...
    'RowNames', {'From Analyitical PSD','From LYAP','From var.m'});

% Pitch damper on:
VAR = [var_on; lyap_var_on; V_on];
SIGMA_on = table(VAR(:, 1), VAR(:, 2), VAR(:, 3), VAR(:, 4), VAR(:, 5),...
    'VariableNames', {'sigma_V','sigma_alpha','sigma_theta','sigma_q','sigma_Nz'},...
    'RowNames', {'From Analyitical PSD','From LYAP','From var.m'});


display(SIGMA_off)
display(SIGMA_on)
%%
error13 = 100*(var_on-V_on)./var_on;

error23 = 100*(lyap_var_on-V_on)./lyap_var_on

error12 = 100*(var_on-lyap_var_on)./var_on

%%
if visualize
    disp('  Press any button to close all windows.')
    pause
end
close all;