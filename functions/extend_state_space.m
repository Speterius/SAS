function ss_ext = extend_state_space(sys, V, c)

% INPUT: STATE SPACE SYSTEM 7x7: [u alpha theta q u_g alpha_g alpha_g_dot]

% OUTPUT: Extended state space with N_z = V/g*(q - alpha_dot) added to C
% and D matrices

    if ~isa(sys, 'ss')
        error('Input error: input has to be a state space object.')
    end

% Load factor rows: N_z = V/g*(q - alpha_dot)
    g = 9.80665;
    N_z_C = V/g * (V/c*[0 0 0 1 0 0 0] - sys.A(2,:));
    N_z_D = V/g * (-sys.B(2,:));
    
% New A, B, C and D matrices:
    A = sys.A;
    B = sys.B;
    C = [sys.C; N_z_C];
    D = [sys.D; N_z_D];
    
% Add to the new output state space:
    ss_ext = ss(A, B, C, D);
    ss_ext.StateName    = {'u','\alpha','\theta','qc/V','u_g','alpha_g','alpha_g_dot'};
    ss_ext.InputName    = {'\delta_{el}','w1','w3'};
    ss_ext.OutputName   = {'u','\alpha','\theta','qc/V','u_g','alpha_g','alpha_g_dot', 'load_factor'};

end