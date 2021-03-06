function [omega, S_V, S_alpha, S_theta, S_q, S_N] = experi_psd(sys, dt, T, seed)

% Based on exampl83.m

% INPUT:    sys:        state space system                      [ss]
%           w:          frequency logspace                      [log rad/s]
%           dt:         time setting deltat                     [s]
%           T_max:      time setting max T                      [s]
%           seed:       RNG seed integer                        [-]

% OUTPUT:   PSD as a function of w freq logspace
%           PSD_pwelch_out - 5 column matrix containing PSD pwelch data
%           w_pwelch_out   - frequency data from the pwelch routine

% Check for input errors:
    if ~isa(sys, 'ss')
        error('Input Error: sys input should be state space object.')
    end
    
% TIME DEFINITION
    t  = 0:dt:T;                % T_max seconds of data
    N = length(t);                  % number of time samples
    
% INPUT VECTOR DEFINITION:
    rng(seed)

    nn = zeros(1,N);                % zero input elevator
    w1 = zeros(1,N);                % zero horizontal gust
    w3 = randn(1,N)/sqrt(dt);       % white noise vertical gust
    
    u  = [nn' w1' w3'];
    
% TIME OUTPUT:
    y = lsim(sys,u,t);
    
% State outputs:
    V = y(:,1);
    alpha = y(:,2);
    theta = y(:,3);
    qc_V = y(:,4);
    N_z = y(:,8);

% COMPUTE PERIODOGRAM AND ESTIMATE PSD
% PERIODOGRAM
    V_f       = dt*fft(V);
    alpha_f   = dt*fft(alpha);
    theta_f   = dt*fft(theta);
    qc_V_f    = dt*fft(qc_V);
    N_z_f     = dt*fft(N_z);

% PSD ESTIMATE
    S_V     = (1/T)*(V_f.*conj(V_f));
    S_alpha = (1/T)*(alpha_f.*conj(alpha_f));
    S_theta = (1/T)*(theta_f.*conj(theta_f));
    S_q     = (1/T)*(qc_V_f.*conj(qc_V_f));
    S_N     = (1/T)*(N_z_f.*conj(N_z_f));

% DEFINE FREQUENCY VECTOR
    fs = 1/dt;     % sample frequency
    omega = 2*pi*fs*(0:(N/2)-1)/N;

% Only positive side:
    S_V         = S_V(1:N/2);
    S_alpha     = S_alpha(1:N/2);
    S_theta     = S_theta(1:N/2);
    S_q         = S_q(1:N/2);
    S_N         = S_N(1:N/2);

end