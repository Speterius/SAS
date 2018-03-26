clear all;
clc;
%%
dt = 0.1;
t = 0:dt:500;
N = size(t, 2);
x = cos(pi/4*t) + .5*randn(size(t));
y = cos(6*pi/4*t) + .5*randn(size(t));
A = [x', y'];

%%

% default
[pAA_def, w_def] = pwelch(A);

% mine
window = 3*N/5;
noverlap = window/2;
fs = 1/dt;
nfft = 50*fs;

[pAA, f] = pwelch(A, window, noverlap, nfft, fs);

% collect
pxx_def = pAA_def(:,1);
pyy_def = pAA_def(:,2);
pxx = pAA(:,1);
pyy = pAA(:,2);

pxx_def_db = 10*log(pxx_def);
pxx_db = 10*log(pxx);

% plot
figure;
plot(w_def/pi, pxx_def_db, '--');
hold on
plot(f, pxx_db);
hold off;
legend('pxx def', 'my pxx')
grid;

disp(size(w_def))
disp(size(w))