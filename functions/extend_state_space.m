function ss_ext = extend_state_space(sys)

% INPUT: STATE SPACE SYSTEM 7x7: [u alpha theta q u_g alpha_g alpha_g_dot]

% OUTPUT: Extended state space with gamma_dot = theta_dot - alpha_dot added to C
% and D matrices

    if ~isa(sys, 'ss')
        error('Input error: input has to be a state space object.')
    end

% alpha and gamma are 2nd and 3rd respectively:
    gamma_dot_C = sys.A(3,:)-sys.A(2,:);
    gamma_dot_D = sys.B(3,:)-sys.B(2,:);
    
% New A, B, C and D matrices:
    A = sys.A;
    B = sys.B;
    C = [sys.C; gamma_dot_C];
    D = [sys.D; gamma_dot_D];

% Add to the new output state space:
    ss_ext = ss(A, B, C, D);
    ss_ext.StateName    = {'u','\alpha','\theta','qc/V','u_g','alpha_g','alpha_g_dot'};
    ss_ext.InputName    = {'\delta_{el}','w1','w3'};
    ss_ext.OutputName   = {'u','\alpha','\theta','qc/V','u_g','alpha_g','alpha_g_dot', 'gamma_dot'};

end