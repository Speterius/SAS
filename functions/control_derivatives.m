%% Cessna Ce500 Citation I
% Stability and Control derivative data:
%LANDING CONFIGURATION
%% 
x_cg            = 0.3;      % [chord length] 
W               = 44675;    % [N]
m               = 4556;     % [kg]
S               = 24.2;     % [m2]
c               = 2.022;    % [m]
b               = 13.36;    % [m]

V               = 51.4;     % [m/s]  
h               = 0;        % [m]
rho             = 1.225;    % [kg/m3]  
mu_c            = 76;       % [-]

mu_b            = 11;       % [-]
K_X2            = 0.012;    % [-]
K_Z2            = 0.037;    % [-]
K_XZ            = 0.002;    % [-]
K_Y2            = 0.980;    % [-]

%% Symetrical
C_X_0           =  0.0;      % [-]
C_X_u           = -0.2173;   % [-]
C_X_alpha       =  0.4692;   % [-]
C_X_q           =  0.0;      % [-]
C_X_delta       =  0.0;      % [-]

C_Z_0           = -1.136;    % [-]
C_Z_u           = -2.272;    % [-]
C_Z_alpha       = -5.13;     % [-]
C_Z_alpha_dot   = -1.405;    % [-]
C_Z_q           = -3.84;     % [-]
C_Z_delta       = -0.6238;   % [-]

C_Z_alpha_g_dot = C_Z_alpha_dot - C_Z_q;

C_m_u           =  0.0;      % [-]
C_m_alpha       = -0.4;      % [-]
C_m_alpha_dot   = -3.615;    % [-]
C_m_q           = -7.35;     % [-]
C_m_delta       = -1.553;    % [-]

C_m_alpha_g_dot = C_m_alpha_dot - C_m_q;

%% Asymetrical
C_Y_beta    = -0.9896;       % [-]
C_Y_p       = -0.087;        % [-]
C_Y_r       =  0.43;         % [-]
C_Y_da      =  0;            % [-]
C_Y_dr      =  0.3037;       % [-]

C_l_beta    = -0.0772;       % [-]
C_l_p       = -0.3415;       % [-]
C_l_r       =  0.2830;       % [-]
C_l_da      = -0.2349;       % [-]
C_l_dr      =  0.0286;       % [-]

C_n_beta    =  0.1628;       % [-]
C_n_p       = -0.0108;       % [-]
C_n_r       = -0.1930;       % [-]
C_n_da      =  0.0286;       % [-]
C_n_dr      = -0.1261;       % [-]

%%
save('cessna_data.mat')

