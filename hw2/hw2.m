clear;clc;close all





sig_A = 3;

sig_B = 8;


1/(1/sig_A^2 - 1/sig_B^2)

sig_A^2 * sig_B^2 / (sig_B^2 - sig_A^2)



%% 5
clear;clc;close all


alpha = .001;

n = 10; % samples

n2 = 2*n;

k_q0 = chi2inv(1 - alpha, n2) / 2;

% q0 = k_q0 *

Pmd = chi2cdf(2*k_q0/4, n2);


fprintf('5: \n\tk_q0: %f\n\tPmd: %f\n', k_q0, Pmd)


%% 2-11
clear;clc;close all

k = 1; % number of measurements

sig = 10^2; % m
d = 10^5 * ones(1, k); % m
y = 10^3; % m

% estimation
J = -1/sig^2 * sum(d.^2 ./ (d.^2 + y^2) - 1);

Jinv = 1/J;

sig_yhat = sqrt(Jinv); % m

n = 10000;
yhats = nan(1, n);
for i = 1:n
    zs = sqrt(d(1)^2 + y^2) + sig*randn(1, k); % measurments, m

    zbar = max(mean(zs), d(1)); % truncate values that cause problems...

    yhat = sqrt(zbar^2 - d(1)^2); % MLE estimate, m

    yhats(i) = yhat;
end % for

histfit(yhats) % plot histogram and fitted normal pdf
pd = fitdist(yhats', 'Normal'); % return fitted normal pdf


fprintf('2-11: \n\tTrue Altitude: %f (m)\n\tSTD of CRLB: %f (m)\n\tObserved Estimator STD: %f (m)\n\tObserved Estimator Mean: %f (m)\n', y, sig_yhat, pd.sigma, pd.mu)













