%% Compare Optimal Neyman-Pearson hypothesis test with hueristics derived test
% Spencer Freeman, 10/21/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 7 of problem set 3 which is highy related to 
% number 1-9 (Bar Shalom) of problem set 1.
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW3-P7_midterm')



%% a

alpha = .01;
P = [1 .5; ...
    .5 2];

Pinv = inv(P);

e = [1; 1];

Pe = Pinv * e;

sig_beta = sqrt(Pe' * P * Pe); % variance of beta
mu_beta  = 0;                  % mean of beta

beta0 = -norminv(alpha/2, mu_beta, sig_beta); % threshold value


% create sample measurements and assess the test
thetas =  -10:.01:10;
for i = 1:length(thetas)
    theta = thetas(i);         % signal
    m = 100;%100e3;                 % number of samples
    w = mvnrnd([0; 0], P, m)'; % random draw noise terms
    z = theta * e + w;         % noisy samples
    b = z' * Pinv * e;         % test statistic for each sample

    accept_H1 = abs(b) >= beta0; % test hypothesis

    pw_beta(i) = sum(accept_H1) / m; % detection rate (power)
    Power_beta(i) = ...
        normcdf(-beta0,   (theta * e)' * Pinv * e, sig_beta) + ...
        1-normcdf( beta0, (theta * e)' * Pinv * e, sig_beta);

end



%% b

bs = linspace(-5, 10, 500); % beta's to evaluate

sig_beta = sqrt(Pe' * P * Pe); % variance of beta
mu_beta  = 0;                  % mean of beta

y0 = normpdf(bs, mu_beta, sig_beta);

theta1 = 4;
mu_beta  = theta1*e'*Pe;                  % mean of beta

y1 = normpdf(bs, mu_beta, sig_beta);



%% plotting 
close all

% Be sure to hand in your acquisition test statistic's formula, 
% its threshold value, and its probability density functions, 
% all with numerical values included where appropriate.

% CDF's of beta and eta
h = figure;
h.WindowStyle = 'Docked';
plot(thetas, pw_beta, 'o', 'Color', "#0072BD"); hold on
plot(thetas, Power_beta, 'LineWidth', 1.5, 'Color', "#D95319")
grid on
title('Part a')
ylabel('Power')
xlabel('Theta')
legend('Observed-Beta', 'Theory-Beta')

% PDF's for beta
h = figure;
h.WindowStyle = 'Docked';
plot(bs, y0, 'LineWidth', 1.5)
hold on
plot(bs, y1, 'LineWidth', 1.5)
grid on
title('Part b')

xline(beta0, 'LineWidth', 1.5)
ylabel('Probability Density')
xlabel('\beta')
legend('Theta = 0', 'Theta = 4', 'Threshold \beta')

fprintf('\n\tThreshold Beta0: %f\n\t1-Sigma Beta: %f\n', beta0, sig_beta)






