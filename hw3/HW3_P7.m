%% Compare Optimal Neyman-Pearson hypothesis test with hueristics derived test
% Spencer Freeman, 10/19/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 7 of problem set 3 which is highy related to 
% number 1-9 (Bar Shalom) of problem set 1.
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW3-P7')



%% a

alpha = .01;
P = [1 .5; .5 1];

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



%% c
% for derivation see HW1_P1-9

k = [1; -.3];

sig_eta = sqrt(k' * P * k); % variance of eta
mu_eta  = 0;                % mean of eta (theta=0)

eta0 = -norminv(alpha/2, mu_eta, sig_eta); % threshold value


%% d

% create sample measurements and assess the test
thetas =  -10:.01:10;
for i = 1:length(thetas)
    theta = thetas(i);         % signal
    m = 100;%100e3;                 % number of samples
    w = mvnrnd([0; 0], P, m)'; % random draw noise terms
    z = theta * e + w;         % noisy samples
    n = k' * z;         % test statistic for each sample

    accept_H1 = abs(n) >= eta0; % test hypothesis

    pw_eta(i) = sum(accept_H1) / m; % detection rate (power)
    Power_eta(i) = ...
        normcdf(-eta0,   (theta * e)' * k, sig_eta) + ...
        1-normcdf( eta0, (theta * e)' * k, sig_eta);

end

%% e

ns = linspace(-3.5, 6.5, 500); % eta's to evaluate

sig_eta = sqrt(k' * P * k); % variance of beta
mu_eta  = 0;                  % mean of beta

y0_eta = normpdf(ns, mu_eta, sig_eta);

theta1 = 4;
mu_eta  = theta1*e'*k;                  % mean of beta

y1_eta = normpdf(ns, mu_eta, sig_eta);


%% f

% The Neyman-Pearson test is more powerful in comparison to heuristic tests.
% This is apparent in the Power coplot where the heuristic test is lower
% than the optimal test around theta=0. The PDF plot of the hueristic also
% shows the need to have a high threshold to reject false positives (the
% PDF's overlap.

%% plotting 
close all

% CDF's of beta and eta
h = figure;
h.WindowStyle = 'Docked';
plot(thetas, pw_beta, 'o', 'Color', "#0072BD"); hold on
plot(thetas, Power_beta, 'LineWidth', 1.5, 'Color', "#D95319")
plot(thetas, pw_eta,  "*", 'Color', "#0072BD")
plot(thetas, Power_eta, '--', 'LineWidth', 1.5, 'Color', "#D95319")
grid on
title('Part a + d')
ylabel('Power')
xlabel('Theta')
legend('Observed-Beta', 'Theory-Beta', 'Observed-Eta', 'Theory-Eta')

% PDF's for beta
h = figure;
h.WindowStyle = 'Docked';
plot(bs, y0, bs, y1)%, thetas, Power)
grid on
title('Part b')

xline(beta0)
ylabel('Probability Density')
xlabel('\beta')
legend('Theta = 0', 'Theta = 4', 'Threshold \beta')

% PDF's for eta
h = figure;
h.WindowStyle = 'Docked';
plot(ns, y0_eta, ns, y1_eta)%, thetas, Power)
grid on
title('Part e')

xline(eta0)
ylabel('Probability Density')
xlabel('\eta')
legend('Theta = 0', 'Theta = 4', 'Threshold \eta')










