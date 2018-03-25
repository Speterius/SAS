n = 0:319;
x = cos(pi/4*n) + randn(size(n));
y = cos(6*pi/4*n) + randn(size(n));
%%
A = [x', y'];

[pAA, w] = pwelch(A);

%%
pxx = pAA(:,1);
pyy = pAA(:,2);
figure;
loglog(w, pxx, w, pyy); 
legend('pxx', 'pyy')