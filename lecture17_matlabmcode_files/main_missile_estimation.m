clear;clc;close all




N_mc = 10;

v0 = 1000; % m/s
el0 = 45*pi/180; % rad
xtrue = [0; v0*cos(el0); 0; v0*sin(el0)];
t1 = 0; % s
delt = 1/10; % s
k = 100; % samples
lradar = 5e3; % m
sigmarho = 1; % m
sigmatheta = 1e-3; % rad




[xerr_mean,P_GN,P_CRLB,Ptypapprox] = ...
    mcsim_gnestmissle(N_mc,xtrue,t1,delt,k,lradar,...
    sigmarho,sigmatheta);


% run 1 MC

xerr_mean = zeros(4,1);
P_GN = zeros(4,4);
J_Fisher = zeros(4,4);

[thist,zhist] = truthmodelmissle(xtrue,t1,delt,k,lradar,...
    sigmarho,sigmatheta);
[xest,Jopt,Ptypapprox,itermflag] = ...
    gnestmissle(thist,zhist,lradar,...
    sigmarho,sigmatheta,0);
if itermflag == 2
    disp(['Warning, the nonlinear least-squares algorithm',...
        ' failed to terminate on iteration ',int2str(ll),'.'])
end
xtil = xtrue - xest;
xerr_mean = xerr_mean + xtil;
P_GN = P_GN + xtil*(xtil');
[J,delxgn,dJdalpha,d2Jdalpha2,P,dJdx] = ...
    jdxgnmissle(xtrue,thist,zhist,lradar,...
    sigmarho,sigmatheta,1);
J_Fisher = J_Fisher + dJdx*(dJdx');


















