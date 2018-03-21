function [A, B, sigmaug_V, sigmaag, V, c] = generate_ac_model(aircraft_data, Lg, sigma)
    %   DESCR:      based on cit2s.m of AE4304 sample codes. 
    %   INPUTS:     aircraft_data:  name of .mat file with stability 
    %                               and control derivatives
    %               Lg:             turbulence scale length
    %               sigma:          turbulence intensity
    %   OUTPUTS:    A and B

    load(aircraft_data)
    sigma_wg = sigma;

    % Turbulence parameters:
    sigmaug_V = sigma_wg/V;
    sigmaag   = sigma_wg/V;

    % Symmetric stability derivatives based on Table 7.1 of reader:

    x_u         = V/c*C_X_u     /2/mu_c;
    x_alpha     = V/c*C_X_alpha /2/mu_c;
    x_theta     = V/c*C_Z_0     /2/mu_c;
    x_q         = V/c*C_X_q     /2/mu_c;
    x_delta     = V/c*C_X_delta /2/mu_c;
    x_ug        = x_u;
    x_ug_d      = 0;
    x_alpha_g   = x_alpha;
    x_alpha_g_d = 0;

    z_u         = V/c*  C_Z_u         /(2*mu_c-C_Z_alpha_dot);
    z_alpha     = V/c*  C_Z_alpha     /(2*mu_c-C_Z_alpha_dot);
    z_theta     = V/c*  (-C_X_0)      /(2*mu_c-C_Z_alpha_dot);
    z_q         = V/c*  (2*mu_c+C_Z_q)/(2*mu_c-C_Z_alpha_dot);
    z_delta     = V/c*  C_Z_delta     /(2*mu_c-C_Z_alpha_dot);
    z_ug        = z_u;
    z_ug_d      = 0;
    z_alpha_g   = z_alpha;
    z_alpha_g_d = V/c*(C_Z_alpha_g_dot)/(2*mu_c-C_Z_alpha_dot);

    m_u         = V/c*  ((C_m_u+C_Z_u*C_m_alpha_dot)         /(2*mu_c-C_Z_alpha_dot))/(2*mu_c*K_Y2);
    m_alpha     = V/c*  ((C_m_alpha+C_Z_alpha*C_m_alpha_dot) /(2*mu_c-C_Z_alpha_dot))/(2*mu_c*K_Y2);
    m_theta     = V/c*  ((-C_X_0*C_m_alpha_dot)              /(2*mu_c-C_Z_alpha_dot))/(2*mu_c*K_Y2);
    m_q         = V/c*  ((C_m_q+C_m_alpha_dot*(2*mu_c-C_Z_q))/(2*mu_c-C_Z_alpha_dot))/(2*mu_c*K_Y2);
    m_delta     = V/c*  ((C_m_delta+C_Z_delta*C_m_alpha_dot) /(2*mu_c-C_Z_alpha_dot))/(2*mu_c*K_Y2);
    m_ug        = m_u;
    m_ug_d      = 0;
    m_alpha_g   = m_alpha;
    m_alpha_g_d = V/c*  ((C_m_alpha_g_dot+C_Z_alpha_g_dot*C_m_alpha_dot) /(2*mu_c-C_Z_alpha_dot))/(2*mu_c*K_Y2);


    % STATE- AND INPUT MATRICES based on eq: 7.109)
    A=[x_u  x_alpha     x_theta     0       x_ug                            x_alpha_g       0;
       z_u  z_alpha     z_theta     z_q     z_ug-z_alpha_g_d*V/Lg*(c/V)     z_alpha_g       z_alpha_g_d*(c/V);
       0        1       0           V/c     0                               0                   0;
       m_u  m_alpha     m_theta     m_q     m_ug-m_ug_d*V/Lg*(c/V)          m_alpha_g       m_alpha_g_d*(c/V);
       0        0       0           0       -V/Lg                           0                   0;
       0        0       0           0       0                               0                   1;
       0        0       0           0       0                               -(V/Lg)^2       -2*V/Lg];

    B=...
     [x_delta       0                                           0;
      z_delta       z_ug_d*(c/V)*sigmaug_V*sqrt(2*V/Lg)         z_alpha_g_d*(c/V)*sigmaag*sqrt(3*V/Lg);
      0             0                                           0;
      m_delta       m_ug_d*(c/V)*sigmaug_V*sqrt(2*V/Lg)         m_alpha_g_d*(c/V)*sigmaag*sqrt(3*V/Lg);
      0             sigmaug_V*sqrt(2*V/Lg)                      0;
      0             0                                           sigmaag*sqrt(3*V/Lg);
      0             0                                           (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3)];


end

