function S_out = analytic_psd(sys, input_nr, output_nr, w)

% Based on exampl83.m

% INPUT:    sys:        state space system                      [ss]
%           input_nr:   index of input variable in B matrix     [-]
%           output_nr:  index of output variable in C matrix    [-]
%           w;          frequency logspace                      [log rad/s]

% OUTPUT:   PSD as a function of w freq logspace

% Check for input errors:
    if ~isa(sys, 'ss')
        error('Input Error: sys input should be state space object.')
    end

% initialize matrices:
    A = sys.A;
    B = sys.B;
    C = sys.C(output_nr, :);
    D = sys.D(output_nr, :);

% get bode data:
    temp = bode(A,B,C,D,input_nr,w); 
    S_out  = temp.*temp;

end