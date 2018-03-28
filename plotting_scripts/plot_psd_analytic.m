
% Velocity
figure('rend','painters','pos',[10 10 900 600]); 
loglog(w_a, S_V_off_a, w_a, S_V_on_a, 'Linewidth', 2);
title('Analytical Velocity PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{VV}'); xlabel('log w'); grid;

% Alpha
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_a, S_alpha_off_a, w_a, S_alpha_on_a, 'Linewidth', 2);
title('Analytical Angle of attack PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{\alpha\alpha}'); xlabel('log w'); grid;

% Theta
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_a, S_theta_off_a, w_a, S_theta_on_a, 'Linewidth', 2);
title('Analytical Pitch angle PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{\theta\theta}'); xlabel('log w'); grid;

% Pitch rate
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_a, S_q_off_a, w_a, S_q_on_a, 'Linewidth', 2);
title('Analytical Pitch rate PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{qq}'); xlabel('log w'); grid;

% Load factor
figure('rend','painters','pos',[10 10 900 600]);  
loglog(w_a, S_N_off_a, w_a, S_N_on_a, 'Linewidth', 2);
title('Analytical Load factor PSD'); legend('Pitch damper OFF','Pitch damper ON')
ylabel('S_{NN}'); xlabel('log w'); grid;