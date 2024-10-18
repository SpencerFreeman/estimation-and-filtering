function [psihist,y1hist,y2hist] = psiy1y2cart(x,thist)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function gives the current heading angle and y1-y2 position
%  of the tricycle cart.
%  
%
%  Inputs:
%
%    x        = [psi0;y10;y20;psidot;vrear], the vector of initial
%             conditions and rates of the tricycle cart at time 0.  psi0
%             is the initial heading angle in rad, y10 is the initial 
%             east position in meters, y20 is the initial north position
%             in meters, psidot is the turn rate in rad/sec, and vrear
%             is the speed of the midpoint between the two rear
%             wheels in m/sec.
%
%    thist    The kx1 vector of time points at which the heading and
%             position is to be output, in seconds.
%  
%  Outputs:
%
%    psihist  The kx1 vector that contains the heading angles at the
%             times in thist, in radians.
%
%    y1hist   The kx1 vector that contains the east positions at the
%             times in thist, in meters.
%
%    y2hist   The kx1 vector that contains the north positions at the
%             times in thist, in meters.
%

%
%  Compute the position of the cart at the times in thist.
%
   k = size(thist,1);
   psihist = zeros(k,1);
   y1hist = zeros(k,1);
   y2hist = zeros(k,1);
   for j = 1:k
      tj = thist(j,1);
      q = x(4)*tj*0.5;
      [sa,dsadq,d2sadq2] = safunct(q);
      pfac = x(5)*tj*sa;
      gamma = x(1) + q;
      cosgamma = cos(gamma);
      singamma = sin(gamma);
      y1hist(j,1) = x(2) + pfac*cosgamma;
      y2hist(j,1) = x(3) + pfac*singamma;
      psihist(j,1) = x(1) + 2*q;
   end