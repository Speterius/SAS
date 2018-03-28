
% smoothing_parameter
MA = 1;

% Velocity
figure('rend','painters','pos',[10 10 900 600]); 
loglog(w_p, smooth(S_V_off_p, MA), w_p, smooth(S_V_on_p, MA), 'Linewidth', 1);
title('Welch PSD Velocity PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{VV}'); xlabel('log w'); grid;

% Alpha
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_p, smooth(S_alpha_off_p, MA), w_p, smooth(S_alpha_on_p, MA), 'Linewidth', 1);
title('Welch PSD Angle of attack PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{\alpha\alpha}'); xlabel('log w'); grid;

% Theta
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_p, smooth(S_theta_off_p, MA), w_p, smooth(S_theta_on_p, MA), 'Linewidth', 1);
title('Welch PSD Pitch angle PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{\theta\theta}'); xlabel('log w'); grid;

% Pitch rate
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_p, smooth(S_q_off_p, MA), w_p, smooth(S_q_on_p, MA), 'Linewidth', 1);
title('Welch PSD Pitch rate PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{qq}'); xlabel('log w'); grid;

% Load factor
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_p, smooth(S_N_off_p, MA), w_p, smooth(S_N_on_p, MA), 'Linewidth', 1);
title('Welch PSD Load factor PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{NN}'); xlabel('log w'); grid;