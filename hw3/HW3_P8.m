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
thetas = 2;% -10:.01:10;
for i = 1:length(thetas)
    theta = thetas(i);         % signal
    m = 100;%100e3;                 % number of samples
    w = mvnrnd([0; 0], P, m)'; % random draw noise terms
    z = theta * e + w;         % noisy samples
    b = z' * Pinv * e;         % test statistic for each sample

    r = nan(m, 1);
    theta_hats = nan(m, 1);
    for j = 1:m
        zu = z(:, j);

        theta_hat = 1/(e'*Pe)*zu'*Pe;
        theta_hats(j) = theta_hat;

        r(j) = 1/2*(zu'*Pinv*zu - ...
            (zu - theta_hat * e)' * Pinv * (zu - theta_hat * e));
    end

    accept_H1_t = r.^2 >= beta0^2; % test hypothesis
    % accept_H1_t = b.^2 >= beta0^2; % test hypothesis
    pw_beta_t(i) = sum(accept_H1_t) / m; % detection rate (power)

    accept_H1 = b.^2 >= beta0^2; % test hypothesis
    % accept_H1 = abs(b) >= beta0; % test hypothesis

    pw_beta(i) = sum(accept_H1) / m; % detection rate (power)
    Power_beta(i) = ...
        normcdf(-beta0,   (theta * e)' * Pinv * e, sig_beta) + ...
        1-normcdf( beta0, (theta * e)' * Pinv * e, sig_beta);

end



%% plotting 
close all

figure; plot(r, '-o', 'LineWidth', .75); hold on; plot(b, '-o'); yline(beta0);grid on

figure;histogram(r*1)
figure;histogram(b.^2)
% 
% figure;histfit(theta_hats)
% pd = fitdist(theta_hats,'Normal')

% CDF's of beta and eta
h = figure;
h.WindowStyle = 'Docked';
plot(thetas, pw_beta, 'o', 'Color', "#0072BD"); hold on
plot(thetas, Power_beta, 'LineWidth', 1.5, 'Color', "#D95319")
% plot(thetas, pw_eta,  "*", 'Color', "#0072BD")
% plot(thetas, Power_eta, '--', 'LineWidth', 1.5, 'Color', "#D95319")
plot(thetas, pw_beta_t, '.', 'Color', "r")

grid on
title('Part a + d')
ylabel('Power')
xlabel('Theta')
legend('Observed-Beta', 'Theory-Beta', 'Observed-Eta', 'Theory-Eta')

% % PDF's for beta
% h = figure;
% h.WindowStyle = 'Docked';
% plot(bs, y0, bs, y1)%, thetas, Power)
% grid on
% title('Part b')
% 
% xline(beta0)
% ylabel('Probability Density')
% xlabel('\beta')
% legend('Theta = 0', 'Theta = 4', 'Threshold \beta')
% 
% % PDF's for eta
% h = figure;
% h.WindowStyle = 'Docked';
% plot(ns, y0_eta, ns, y1_eta)%, thetas, Power)
% grid on
% title('Part e')
% 
% xline(eta0)
% ylabel('Probability Density')
% xlabel('\eta')
% legend('Theta = 0', 'Theta = 4', 'Threshold \eta')










