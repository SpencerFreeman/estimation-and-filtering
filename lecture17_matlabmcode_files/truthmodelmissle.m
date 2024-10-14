function [thist,zhist] = truthmodelmissle(xtrue,t1,delt,k,lradar,...
                                            sigmarho,sigmatheta)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function generates radar "measurements" via truth-model
%  simulation of the missle tracking problem.
%  
%
%  Inputs:
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
%    thist       = t1 + [0:(k-1)]'*delt, the kx1 vector of radar
%                measurement times.
%
%    zhist       = [rho_1,theta_1;rho_2,theta_2;rho_3,theta_3;...;...
%                rho_k,theta_k], the kx2 time history of the simulated
%                radar measurements.
%

%
%  Set up the output arrays.
%
   thist = t1 + [0:(k-1)]'*delt;
   zhist = zeros(k,2);
%
%  Loop through the measurement times and compute the "perfect"
%  radar measurements, i.e., the radar measurements without
%  noise.
%
   for j = 1:k
      hj = hjmissle(xtrue,thist(j,1),lradar,0,0);
      zhist(j,:) = hj';
   end
%
%  Add measurement noise to the outputs.
%
   zhist(:,1) =   zhist(:,1) + randn(k,1)*sigmarho;
   zhist(:,2) =   zhist(:,2) + randn(k,1)*sigmatheta;