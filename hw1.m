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










