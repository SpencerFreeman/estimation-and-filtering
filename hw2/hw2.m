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









