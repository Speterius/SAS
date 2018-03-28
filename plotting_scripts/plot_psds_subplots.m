%% WITHOUT SMOOTHING
% V) subplot with pitch damper on/off
MA = 1;

figure('rend','painters','pos',[10 10 900 600]); 

subplot(2,1,1);
loglog(w_e, S_V_off_e, '-.'); hold on;
loglog(w_p, S_V_off_p, 'Linewidth', 2); 
loglog(w_a, S_V_off_a, 'Linewidth', 1);
title('S_{VV}(\omega) with Pitch damper off'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{VV}'); xlabel('log w'); grid;

subplot(2,1,2)
loglog(w_e, S_V_on_e, '-.'); hold on;
loglog(w_p, S_V_on_p, 'Linewidth', 2);
loglog(w_a, S_V_on_a, 'Linewidth', 1)
title('S_{VV}(\omega) with Pitch damper on'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{VV}'); xlabel('log w'); grid;

%% ALPHA) subplots

figure('rend','painters','pos',[10 10 900 600]); 

subplot(2,1,1);
loglog(w_e, S_alpha_off_e, '-.'); hold on;
loglog(w_p, S_alpha_off_p, 'Linewidth', 2); 
loglog(w_a, S_alpha_off_a, 'Linewidth', 1);
title('S_{\alpha\alpha}(\omega) with Pitch damper OFF'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{\alpha\alpha}'); xlabel('log w'); grid;

subplot(2,1,2)
loglog(w_e, S_alpha_on_e, '-.'); hold on;
loglog(w_p, S_alpha_on_p, 'Linewidth', 2);
loglog(w_a, S_alpha_on_a, 'Linewidth', 1)
title('S_{\alpha\alpha}(\omega) with Pitch damper ON'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{\alpha\alpha}'); xlabel('log w'); grid;

%% THETA) subplots

figure('rend','painters','pos',[10 10 900 600]); 

subplot(2,1,1);
loglog(w_e, S_theta_off_e, '-.'); hold on;
loglog(w_p, S_theta_off_p, 'Linewidth', 2); 
loglog(w_a, S_theta_off_a, 'Linewidth', 1);
title('S_{\theta\theta}(\omega) with Pitch damper OFF'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{\theta\theta}'); xlabel('log w'); grid;

subplot(2,1,2)
loglog(w_e, S_theta_on_e, '-.'); hold on;
loglog(w_p, S_theta_on_p, 'Linewidth', 2);
loglog(w_a, S_theta_on_a, 'Linewidth', 1)
title('S_{\theta\theta}(\omega) with Pitch damper ON'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{\theta\theta}'); xlabel('log w'); grid;

%% q) subplots

figure('rend','painters','pos',[10 10 900 600]); 

subplot(2,1,1);
loglog(w_e, S_q_off_e, '-.'); hold on;
loglog(w_p, S_q_off_p, 'Linewidth', 2); 
loglog(w_a, S_q_off_a, 'Linewidth', 1);
title('S_{qq}(\omega) with Pitch damper OFF'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{qq}'); xlabel('log w'); grid;

subplot(2,1,2)
loglog(w_e, S_q_on_e, '-.'); hold on;
loglog(w_p, S_q_on_p, 'Linewidth', 2);
loglog(w_a, S_q_on_a, 'Linewidth', 1)
title('S_{qq}(\omega) with Pitch damper ON'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{qq}'); xlabel('log w'); grid;

%% Nz) subplots

figure('rend','painters','pos',[10 10 900 600]); 

subplot(2,1,1);
loglog(w_e, S_N_off_e, '-.'); hold on;
loglog(w_p, S_N_off_p, 'Linewidth', 2); 
loglog(w_a, S_N_off_a, 'Linewidth', 1);
title('S_{NN}(\omega) with Pitch damper OFF'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{NN}'); xlabel('log w'); grid;

subplot(2,1,2)
loglog(w_e, S_N_on_e, '-.'); hold on;
loglog(w_p, S_N_on_p, 'Linewidth', 2);
loglog(w_a, S_N_on_a, 'Linewidth', 1)
title('S_{NN}(\omega) with Pitch damper ON'); legend({'FFT Experimental','PWELCH Experimental','Analytical'},'FontSize',12)
ylabel('S_{NN}'); xlabel('log w'); grid;