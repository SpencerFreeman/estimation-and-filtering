function [xerr_mean,P_NEWTON,P_CRLB,Ptypapprox] = ...
                     mcsim_newtonestmissle(N_mc,xtrue,t1,delt,k,lradar,...
                                           sigmarho,sigmatheta)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function does Monte Carlo simulation of the missle 
%  estimation problem using the Newton estimator
%  
%
%  Inputs:
%
%    N_mc        The number of Monte-Carlo simulations to do.
%
%    xtrue       = [y10true;v10true;y20true;v20true], the true vector of 
%                initial conditions of the missle at time 0.  y10true is
%                the true initial horizontal position in meters, v10true
%                is the true initial horizontal velocity in m/sec, y20true 
%                is the true initial vertical position in meters, and 
%                v20true is the true initial vertical velocity in m/sec.
%
%    t1          The time in seconds of the first radar measurement.
%
%    delt        The time in seconds between radar samples.
%
%    k           The total number of radar samples.
%
%    lradar      The horizontal position of the radar station, in meters.
%                The radar station is assumed to be located at zero
%                vertical position.
%
%    sigmarho    The measurement error standard deviation of the
%                radar range measurement, in meters.
%
%    sigmatheta  The measurement error standard deviation of the
%                radar elevation measurement, in radians.
%  
%  Outputs:
%
%    xerr_mean   The 4x1 mean error vector from the N_mc Monte-Carlo
%                simulations.
%
%    P_NEWTON    The 4x4 Mean-Square estimation error matrix (includes
%                any non-zero mean error).
%
%    P_CRLB      The Cramer-Rao lower bound for the estimation error
%                covariance matrix (this is an approximation
%                that approaches the exact Cramer-Rao lower bound
%                as N_mc approaches infinity.
%
%    Ptypapprox  A typical Newton approximation of the
%                estimation error covariance matrix.
%

%
%  Set up the first two outputs for storing interim results.
%
   xerr_mean = zeros(4,1);
   P_NEWTON = zeros(4,4);
   J_Fisher = zeros(4,4);
%
%  Run the Monte Carlo Simulations.
%
   for ll = 1:N_mc
      [thist,zhist] = truthmodelmissle(xtrue,t1,delt,k,lradar,...
                                            sigmarho,sigmatheta);
      [xest,Jopt,Ptypapprox,itermflag] = ...
                     newtonestmissle(thist,zhist,lradar,...
                                     sigmarho,sigmatheta,0);
      if itermflag >= 2
         disp(['Warning, the nonlinear least-squares algorithm',...
               ' failed to terminate well on iteration ',int2str(ll),'.'])
      end
      xtil = xtrue - xest;
      xerr_mean = xerr_mean + xtil;
      P_NEWTON = P_NEWTON + xtil*(xtil');
      [J,delxnewton,dJdalpha,d2Jdalpha2,P,dJdx,d2Jdx2,iflagHess] = ...
                     jdxnewtonmissle(xtrue,thist,zhist,lradar,...
                                     sigmarho,sigmatheta,1);
      J_Fisher = J_Fisher + d2Jdx2;
   end
%
%  Average the outputs of the iterations.
%
   factorN_mc = 1/N_mc;
   xerr_mean = xerr_mean*factorN_mc;
   P_NEWTON = P_NEWTON*factorN_mc;
   J_Fisher = J_Fisher*factorN_mc;
   P_CRLB = inv(J_Fisher);