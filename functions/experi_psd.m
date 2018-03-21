function [omega, S_V, S_alpha, S_theta, S_q, S_N] = experi_psd(sys, dt, T, seed)

% Based on exampl83.m

% INPUT:    sys:        state space system                      [ss]
%           w:          frequency logspace                      [log rad/s]
%           dt:         time setting deltat                     [s]
%           T_max:      time setting max T                      [s]
%           seed:       RNG seed integer                        [-]

% OUTPUT:   PSD as a function of w freq logspace

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
    V       = dt*fft(V);
    alpha   = dt*fft(alpha);
    theta   = dt*fft(theta);
    qc_V    = dt*fft(qc_V);
    N_z     = dt*fft(N_z);

% PSD ESTIMATE
    S_V     = (1/T)*( V.*conj(V));
    S_alpha = (1/T)*(  alpha.*conj(alpha));
    S_theta = (1/T)*(    theta.*conj(theta));
    S_q     = (1/T)*(    qc_V.*conj(qc_V));
    S_N     = (1/T)*(N_z.*conj(N_z));

% DEFINE FREQUENCY VECTOR
    fs = 1/dt;     % sample frequency
    omega = 2*pi*fs*(0:(N/2)-1)/N;

% CUT PSDS IN HALF:
    S_V         = S_V(1:N/2);
    S_alpha     = S_alpha(1:N/2);
    S_theta     = S_theta(1:N/2);
    S_q         = S_q(1:N/2);
    S_N         = S_N(1:N/2);
    
end