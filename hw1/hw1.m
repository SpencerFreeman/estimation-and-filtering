clear;clc;close all





%% 1


S = randi(10, 2, 2);

% S = [ 1,     6; ...
%      3    10];


S+S'


P = (S + S')/2;

x = randi(10, 2, 1);

x'*S*x
x'*P*x


%% 2



%% 1-7
clear;clc;close all

n = 10000;

nx = 2;
ny = 5;

nz = 3;

xbar = randi(10, nx, 1);
ybar = randi(10, ny, 1);

sx = randi(10, nx, 1);
sy = randi(10, ny, 1);

x = randn(nx, n).*sx + xbar;
y = randn(ny, n).*sy + ybar;

A = randi(10, nz, nx);
B = randi(10, nz, ny);
C = randi(10, nz, 1);

z = A*x + B*y + C;

% 1
zbar = mean(z')'
zbar_true = A*xbar + B*ybar + C

P = cov([x; y]');

Pxx = P(1:nx,       1:nx);
Pyy = P(nx + 1:end, nx + 1:end);
Pxy = P(1:nx,       nx + 1:end);
Pyx = Pxy';

% 2
Pz = cov(z')
Pz_true = A*Pxx*A' + A*Pxy*B' + B*Pyx*A' + B*Pyy*B'


%% 1-9
clear;clc;close all

alpha = .01;
P = [1 .5; .5 1];

Pinv = inv(P);

e = [1; 1];

Pe = Pinv * e;

sig_beta = sqrt(Pe' * P * Pe); % variance of beta
mu_beta  = 0;                  % mean of beta

beta0 = -norminv(alpha/2, mu_beta, sig_beta); % threshold value


% create sample measurements and assess the test
thetas = linspace(-5, 5, 100);
for i = 1:length(thetas)
    theta = thetas(i);         % signal
    m = 100;%100e3;                 % number of samples
    w = mvnrnd([0; 0], P, m)'; % random draw noise terms
    z = theta * e + w;         % noisy samples
    b = z' * Pinv * e;         % test statistic for each sample

    accept_H1 = abs(b) >= beta0; % test hypothesis

    pw(i) = sum(accept_H1) / m; % detection rate (power)
    Power(i) = ...
        normcdf(-beta0,   (theta * e)' * Pinv * e, sig_beta) + ...
        1-normcdf( beta0, (theta * e)' * Pinv * e, sig_beta);

end


figure
plot(thetas, pw, 'o', thetas, Power)
grid on
legend('Observed', 'Theory')




%%
x = linspace(-5, 5, 100);
y = normpdf(x, mu_beta, sig_beta);

figure('WindowStyle', 'Docked')
plot(x, y)
grid on
hold on
xline(beta0)
xline(-beta0)



figure;plot(sort(accept_H1));grid on
figure;plot(b);grid on
figure;histogram(b)















