clear;clc;close all



k = 1e6;

xbar = 10;
Pxx = 2*1;
Pww = 7;

w = randn(1, k)*sqrt(Pww);

x = xbar + randn(1, k)*sqrt(Pxx);

z = w + x;
y = z.^2;

ybar = Pxx + Pww + xbar^2;
ybar_obs = mean(y);

Pxy = 3*xbar*Pxx + xbar^3 + xbar*Pww - xbar*ybar;
Pyy = 3*Pxx^2 + 6*xbar^2*Pxx + xbar^4 + 6*(Pxx + xbar^2)*Pww + 3*Pww^2 - ybar^2;
Pyy_obs = var(y);

xhat = xbar + Pxy*inv(Pyy)*(y - ybar);

figure
plot(x, 'o'); hold on; grid on
plot(xhat, 'x')
plot(sqrt(y), '*')

MSE = Pxx - Pxy*inv(Pyy)*Pxy
MSE_obs = mean((xhat - x).^2)
MSE_obs_t = mean((sqrt(y) - x).^2)
MSE_obs_t = mean((xbar - x).^2)


% mean(x.^3)
% 3*xbar*Pxx + xbar^3
% 
% mean(x.^4)
% 3*Pxx^2 + 6*xbar^2*Pxx + xbar^4
% 
% mean(x.*w.^2)
% xbar*Pww
% 
% mean(x.^2.*w.^2)
% (Pxx + xbar^2)*Pww
% 
% mean(w.^4)
% 3*Pww^2

% figure;histogram(x)
% figure;histogram(x.^3)











