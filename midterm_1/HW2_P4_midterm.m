%% Write a Matlab function to solve the weighted least-squares problem
% Spencer Freeman, 10/21/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 4 of problem set 2
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW2-P4_midterm')


% Assignment 2, Number 4, except add 0.25 to each element of z and add 1.2 
% multiplied by the identity matrix to R before solving the problem

z  =  [ -45.1800;...
    1.7900;...
    -31.3800;...
    26.7700;...
    27.6400] + .25;

H  =  [ -4.9300, -1.3100, -1.5900;...
    13.2600, 9.7100, 30.7000;...
    -17.0800, -11.9100, -12.1300;...
    -24.0300, -2.9900, -26.9500;...
    -2.4000, -8.7000, 9.3900];

R  =  [ 5.9700, -0.9200, -1.1800, -7.0600, -1.7900;...
    -0.9200, 3.4500, 1.7100, -0.6000, -4.0500;...
    -1.1800, 1.7100, 1.1900, 0.5600, -1.6700;...
    -7.0600, -0.6000, 0.5600, 9.9200, 4.8500;...
    -1.7900, -4.0500, -1.6700, 4.8500, 6.8700] + 1.2*eye(5);

Ra = chol(R);
Rainv = inv(Ra);
za = Rainv'*z;
Ha = Rainv'*H;
[Qb, Rb0] = qr(Ha);
ind = find(Rb0(:, end) ~= 0, 1, 'last');

zb = Qb'*za;
zbc = zb(1:ind);
Rb = Rb0(1:ind, :);
xhat = inv(Rb)*zbc;

Rinv = inv(R);
sol = norm(-H'*Rinv*(z - H*xhat)) / norm(-H'*Rinv*z);

fprintf('\n\txhat: %1.4f %1.4f %1.4f\n\ttol: %e\n', xhat, sol)
















