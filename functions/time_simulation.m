function [t, V, alpha, theta, qc_V, N_z] = time_simulation(sys, dt, T_max, seed, V_trim, c)

% Based on exampl71.m

% INPUT:    sys:    state space system  [ss]
%           dt:     time_step           [Hz]
%           T_max:  max time            [s]
%           seed:   seed from main      [-]
%           V_trim: trim speed          [m/s]
%           c:      chord length        [m]

% OUTPUT:   t:      time data           [s]
%           + 4 longitudinal states aircraft states
%           + N_z load factor           [g units]

% Check for input errors:
if ~isa(sys, 'ss')
    error('Input Error: sys input should be state space object.')
end
    
% TIME DEFINITION
    t  = 0:dt:T_max;                % T_max seconds of data
    N = length(t);                  % number of time samples
    
% INPUT VECTOR DEFINITION:
    rng(seed)

    nn = zeros(1,N);                % zero input elevator
    w1 = zeros(1,N);                % zero horizontal gust
    w3 = randn(1,N)/sqrt(dt);       % white noise vertical gust
    
    u  = [nn' w1' w3'];

% TIME OUTPUT:
    y = lsim(sys,u,t);
    
% Select relevant states and add load factor + change units to Si
    V = y(:,1)      *V_trim + V_trim;
    alpha = y(:,2)  *180/pi;
    theta = y(:,3)  *180/pi;
    qc_V = y(:,4)   *V_trim*180/pi/c;
    N_z = y(:,8);
end