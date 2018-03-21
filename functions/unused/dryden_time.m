function [t, wg] = dryden_time(sigma, Lg, V, T, dt)
    % DESCR:    based on exampl62.m of AE4304 example matlab codes
    %
    % INPUTS:   sigma:  turbulence intensity    [m/s]
    %           Lg:     turbulence scale length [m]
    %           V:      airspeed                [m/s]
    %           T:      signal length           [s]
    %           dt:     time resolution         [s]
    % OUTPUTS:  t:      timearray
    %           wg:     simulated Dryden gust   [m/s]

    % Define time basis          
    t  = 0:dt:T;
    N = length(t);          

    % White noise input
    w = randn(1,N)/sqrt(dt);
    
    % Forming filter characteristics equation (6.41)
    rat = V/Lg;
    A = [0          1;
        -rat^2      -2*rat];
    B = sigma*[sqrt(3*rat);
        (1-2*sqrt(3))*sqrt((rat^3))];
    C = [1 0];
    D = 0;
    
    % Output turbulence velocity
    wg = lsim(A,B,C,D,w,t);
end