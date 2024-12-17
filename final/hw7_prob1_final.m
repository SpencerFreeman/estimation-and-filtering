%% Fix the code to compute the required partial derivatives.
% Spencer Freeman, 12/17/2024
% AOE 5784, Estimation and Filtering
%
% This script solves number 1 of problem set 7 
% -------------------------------------------------------------------------
clear;clc;close all

disp('HW7-P1-Final')

format long

% Test your code using the supplied test function "fscript_ts01.m". 
% Test the function by numerically integrating from a random initial 
% condition and using a random process noise vector. Integrate over a 
% time span of 3 seconds, and use 120 4th-order Runge-Kutta numerical 
% integration steps. Compare your results with the exact results, which 
% you can computed as outlined in the initial comments section of 
% "fscript_ts01.m". Test the results again, but this time use only 60 
% 4th-order Runge-Kutta steps for the same inputs. Does the error for 
% this second case change as you expect it to change in comparison with 
% the error for the first numerical integration case?

tk = 0; % s
tkp1 = 3; % s

% xk = rand(4, 1);
% vk = rand(3, 1);
% uk = [];
xk = [-0.40; 0.85; -0.60; -1.65]; 
uk = [];
vk = [-0.77; 1.30; 1.65];
idervflag = true;
fscriptname = 'fscript_ts01';

[~, A, D] = fscript_ts01(tk,xk,uk,vk,idervflag);

[dfprinted_dxk_true, dfprinted_dvk_true] = c2d(A, D, (tkp1 - tk));
fprinted_true = dfprinted_dxk_true*xk + dfprinted_dvk_true*vk;

%% many steps
nRK = 120;

[fprinted, dfprinted_dxk, dfprinted_dvk] = ...
    ...
    c2dnonlinear(xk,uk,vk,tk,tkp1,nRK,fscriptname,idervflag);

disp(string(nRK) + " RK Steps")

disp('State Estimate:')
% fprinted_true
fprinted
disp('Partial Derivative wrt xk Estimate:')
% dfprinted_dxk_true
dfprinted_dxk
disp('Partial Derivative wrt vk Estimate:')
% dfprinted_dvk_true
dfprinted_dvk


%% few steps
nRK = 60;

[fprinted, dfprinted_dxk, dfprinted_dvk] = ...
    ...
    c2dnonlinear(xk,uk,vk,tk,tkp1,nRK,fscriptname,idervflag);

disp(string(nRK) + " RK Steps")

disp('State Estimate:')
% fprinted_true
fprinted
disp('Partial Derivative wrt xk Estimate:')
% dfprinted_dxk_true
dfprinted_dxk
disp('Partial Derivative wrt vk Estimate:')
% dfprinted_dvk_true
dfprinted_dvk

disp('Error in state estimate is 10 times larger halving the number of steps. Error in derivatives is nearly doubled.')


























