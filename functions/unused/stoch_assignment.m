


%Get system dynamics
%QUESTION 1 INCLUDED IN cit2a.m
cit2a;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QUESTION 2 FROM HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME AXIS AND INPUT VECTOR DEFINITION
dt = 0.05; T  = 60; t = [0:dt:T]; N = length(t);
nn = zeros(1,N);

% TURBULENCE INPUTS
u_g = randn(1,N)/sqrt(dt);    % sqrt(dt) because of lsim characteristics
w_g = randn(1,N)/sqrt(dt);

% INPUT VECTORS
u1 = [nn' nn' u_g' nn'  nn'];     
u3 = [nn' nn' nn'  w_g'  nn'];

% RESPONSE to u_g
y1 = lsim(sys,u1,t);
% RESPONSE to w_g
y3 = lsim(sys,u3,t);

% RESPONSE TO u_g
clf;
subplot(3,2,1); plot(t,y1(:,1));  
ah=gca; set(ah,'Fontsize',14);
axis([0,60,-0.02,0.02]); xlabel('time (s)'); ylabel('beta');
subplot(3,2,2); plot(t,y1(:,2));  
ah=gca; set(ah,'Fontsize',14);
axis([0,60,-0.2,0.2]); xlabel('time (s)'); ylabel('phi');
subplot(3,2,3); plot(t,y1(:,3)); 
axis([0,60,-0.01,0.01]); xlabel('time (s)'); ylabel('pb/2V');
subplot(3,2,4); plot(t,y1(:,4)); 
axis([0,60,-0.005,0.005]); xlabel('time (s)'); ylabel('rb/2V');
subplot(3,2,5); plot(t,y1(:,5));
axis([0,60,-1,1]); xlabel('time (s)' ); ylabel('ay');
disp("Plotting time response to u_g, press button to continue..");
pause;

% RESPONSE TO w_g
clf;
subplot(3,2,1); plot(t,y3(:,1));
axis([0,60,-0.1,0.1]); xlabel('time, s'); ylabel('beta');
subplot(3,2,2); plot(t,y3(:,2)); 
axis([0,60,-0.1,0.05]);xlabel('time, s'); ylabel('phi');
subplot(3,2,3); plot(t,y3(:,3));
axis([0,60,-0.01,0.01]);xlabel('time, s'); ylabel('pb/2V');
subplot(3,2,4); plot(t,y3(:,4)); 
axis([0,60,-0.01,0.01]); xlabel('time, s'); ylabel('rb/2V');
subplot(3,2,5); plot(t,y3(:,5));
axis([0,60,-1,1]); xlabel('time (s)' ); ylabel('ay');
disp("Plotting time response to w_g, press button to continue..");
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QUESTION 3 FROM HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE FREQUENCY VECTOR
w = logspace(-2,2,300);

% COMPUTE ANALYTIC POWER SPECTRAL DENSITIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOT SURE ABOUT THIS YET
temp = bode(A1,B,C(1,:),D(1,:),3,w); Sbeta  = temp.*temp;
temp = bode(A1,B,C(2,:),D(2,:),3,w); Sphi   = temp.*temp;
temp = bode(A1,B,C(3,:),D(3,:),3,w); Spp    = temp.*temp;
temp = bode(A1,B,C(4,:),D(4,:),3,w); Srr    = temp.*temp;
temp = bode(A1,B,C(5,:),D(5,:),3,w); Say    = temp.*temp;

Sxx1  = [Sbeta Sphi Spp Srr Say];

temp = bode(A1,B,C(1,:),D(1,:),4,w); Sbeta  = temp.*temp;
temp = bode(A1,B,C(2,:),D(2,:),4,w); Sphi   = temp.*temp;
temp = bode(A1,B,C(3,:),D(3,:),4,w); Spp    = temp.*temp;
temp = bode(A1,B,C(4,:),D(4,:),4,w); Srr    = temp.*temp;
temp = bode(A1,B,C(5,:),D(5,:),4,w); Say    = temp.*temp;

Sxx3  = [Sbeta Sphi Spp Srr Say];

% COMPUTE PSDS USING EXPERIMENTAL DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE PERIODOGRAM AND ESTIMATE PSD
% PERIODOGRAM: 1 for gust in u, 3 for gust in w 
temp1   = y1(:,1);
BETA1  = dt*fft(temp1);
temp3   = y3(:,1);
BETA3  = dt*fft(temp3);

temp1    = y1(:,2);
temp3   = y3(:,2);
PHI1   = dt*fft(temp1);
PHI3   = dt*fft(temp3);

temp1    = y1(:,3);
temp3   = y3(:,3);
P1   = dt*fft(temp1);
P3   = dt*fft(temp3);

temp1    = y1(:,4);
temp3   = y3(:,4);
R1   = dt*fft(temp1);
R3   = dt*fft(temp3);

temp1    = y1(:,5);
temp3   = y3(:,5);
AY1   = dt*fft(temp1);
AY3   = dt*fft(temp3);

% PSD ESTIMATE
Pbeta1  = (1/T)*( BETA1.*conj(BETA1));
Pbeta3  = (1/T)*( BETA3.*conj(BETA3));
Pphi1   = (1/T)*(  PHI1.*conj(PHI1));
Pphi3   = (1/T)*(  PHI3.*conj(PHI3));
Pp1     = (1/T)*(    P1.*conj(P1));
Pp3     = (1/T)*(    P3.*conj(P3));
Pr1     = (1/T)*(    R1.*conj(R1));
Pr3     = (1/T)*(    R3.*conj(R3));
Pay1 = (1/T)*(AY1.*conj(AY1));
Pay3 = (1/T)*(AY3.*conj(AY3));

% AVERAGE PSDS USING A SMOOTHING FILTER
Pbeta1s  = smooth(Pbeta1);
Pbeta3s  = smooth(Pbeta3);
Pphi1s   = smooth(Pphi1);
Pphi3s   = smooth(Pphi3);
Pp1s     = smooth(Pp1);
Pp3s     = smooth(Pp3);
Pr1s     = smooth(Pr1);
Pr3s     = smooth(Pr3);
Pay1s = smooth(Pay1);
Pay3s = smooth(Pay3);


% DEFINE FREQUENCY VECTOR
fs = 1/dt;     % sample frequency
omega = 2*pi*fs*(0:(N/2)-1)/N;

% PLOT ANALYTIC AND ESTIMATED PSDS IN ONE PLOT for u_g
subplot(3,2,1); loglog(w,Sxx1(:,1),'--',omega,Pbeta1(1:floor(N/2)));
xlabel('omega [rad/s]'); ylabel('Sbeta');
subplot(3,2,2); loglog(w,Sxx1(:,2),'--',omega,Pphi1(1:floor(N/2)));
axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Sphi');
subplot(3,2,3); loglog(w,Sxx1(:,3),'--',omega,Pp1(1:floor(N/2)));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Spp');
subplot(3,2,4); loglog(w,Sxx1(:,4),'--',omega,Pr1(1:floor(N/2)));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Srr');
subplot(3,2,5); loglog(w,Sxx1(:,5),'--',omega, Pay1(1:floor(N/2)));
axis(10.^[-2 2 -8 2]); xlabel('omega [rad/s]'); ylabel('Say');
disp("Plotting analytic and experimental power spectral density for input u_g. Press any button to continue..");
pause; 

% PLOT ANALYTIC AND SMOOTHED ESTIMATED PSDS IN ONE PLOT for u_g
subplot(3,2,1); loglog(w,Sxx1(:,1),'--',omega,Pbeta1s(1:floor(N/2)));
xlabel('omega [rad/s]'); ylabel('Sbeta');
subplot(3,2,2); loglog(w,Sxx1(:,2),'--',omega,Pphi1s(1:floor(N/2)));
axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Sphi');
subplot(3,2,3); loglog(w,Sxx1(:,3),'--',omega,Pp1s(1:floor(N/2)));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Spp');
subplot(3,2,4); loglog(w,Sxx1(:,4),'--',omega,Pr1s(1:floor(N/2)));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Srr');
subplot(3,2,5); loglog(w,Sxx1(:,5),'--',omega, Pay1s(1:floor(N/2)));
axis(10.^[-2 2 -8 2]); xlabel('omega [rad/s]'); ylabel('Say');
disp("Plotting analytic and smoothed experimental power spectral density for input u_g. Press any button to continue..");
pause; 

% PLOT ANALYTIC AND ESTIMATED PSDS IN ONE PLOT FOR w_g
subplot(3,2,1); loglog(w,Sxx3(:,1),'--',omega,Pbeta3(1:floor(N/2))); 
axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Sbeta');
subplot(3,2,2); loglog(w,Sxx3(:,2),'--',omega,Pphi3(1:floor(N/2)));
axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Sphi');
subplot(3,2,3); loglog(w,Sxx3(:,3),'--',omega,Pp3(1:floor(N/2)));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Spp');
subplot(3,2,4); loglog(w,Sxx3(:,4),'--',omega,Pr3(1:floor(N/2)));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Srr');
subplot(3,2,5); loglog(w,Sxx3(:,5),'--',omega, Pay3(1:floor(N/2)));
axis(10.^[-2 2 -8 2]); xlabel('omega [rad/s]'); ylabel('Say');
disp("Plotting analytic and experimental power spectral density for input w_g. Press any button to continue..");
pause;

% PLOT ANALYTIC AND SMOOTHED ESTIMATED PSDS IN ONE PLOT FOR w_g
subplot(3,2,1); loglog(w,Sxx3(:,1),'--',omega,Pbeta3s(1:floor(N/2))); 
axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Sbeta');
subplot(3,2,2); loglog(w,Sxx3(:,2),'--',omega,Pphi3s(1:floor(N/2)));
axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Sphi');
subplot(3,2,3); loglog(w,Sxx3(:,3),'--',omega,Pp3s(1:floor(N/2)));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Spp');
subplot(3,2,4); loglog(w,Sxx3(:,4),'--',omega,Pr3s(1:floor(N/2)));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Srr');
subplot(3,2,5); loglog(w,Sxx3(:,5),'--',omega, Pay3s(1:floor(N/2)));
axis(10.^[-2 2 -8 2]); xlabel('omega [rad/s]'); ylabel('Say');
disp("Plotting analytic and experimental power spectral density for input w_g. Press any button to continue..");
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QUESTION 4 FROM HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get variance using the analytic PSD (Sxx1 and Sxx3)
%Use is made of the equations on slide 18 of chapter 8 

%Get variance using Lyapunov

%Get variance using time domain data 
V1_tdata = var(y1);
V2_tdata = var(y3);



function s = smooth(n)
    temp = [];
    for i = 2:(length(n)-1)    
        temp = [temp, 0.25*n(i-1) + 0.25 * n(i+1) + 0.5*n(i)];
    end
    temp = [0, temp, n(length(n))];
    s=temp'; 
end
