
% smoothing_parameter
MA = 80;

% Velocity
figure('rend','painters','pos',[10 10 900 600]); 
loglog(w_e, smooth(S_V_off_e, MA), w_e, smooth(S_V_on_e, MA), 'Linewidth', 1);
title('FFT Periodogram Velocity PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{VV}'); xlabel('log w'); grid;

% Alpha
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_e, smooth(S_alpha_off_e, MA), w_e, smooth(S_alpha_on_e, MA), 'Linewidth', 1);
title('FFT Periodogram Angle of attack PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{\alpha\alpha}'); xlabel('log w'); grid;

% Theta
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_e, smooth(S_theta_off_e, MA), w_e, smooth(S_theta_on_e, MA), 'Linewidth', 1);
title('FFT Periodogram Pitch angle PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{\theta\theta}'); xlabel('log w'); grid;

% Pitch rate
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_e, smooth(S_q_off_e, MA), w_e, smooth(S_q_on_e, MA), 'Linewidth', 1);
title('FFT Periodogram Pitch rate PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{qq}'); xlabel('log w'); grid;

% Load factor
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_e, smooth(S_N_off_e, MA), w_e, smooth(S_N_on_e, MA), 'Linewidth', 1);
title('FFT Periodogram Load factor PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{NN}'); xlabel('log w'); grid;