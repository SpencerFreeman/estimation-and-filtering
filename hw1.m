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













