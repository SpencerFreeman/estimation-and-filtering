function xvecguess = fstgessunlat(zhist,lpost)
%
%  Copyright (c) 2022 Mark L. Psiaki.  All rights reserved.  
%
%  "And whatever you do, whether in word or deed, do it all in the name of 
%   the Lord Jesus, giving thanks to God the Father through Him." 
%   (Colossians 3:17).
% 
%  This function derives a crude first guess for the sun-lattitude/
%  user-latitude/user-azimuth problem.  It does this by using a straight-
%  line fit to the data and by guessing that the sun latitude is zero.
%
%  Note: this assumes that the initial azimuth is not more than about
%  45 deg off so that there is no sign uncertainty about north.
%  
%
%  Inputs:
%
%    zhist       = [xshadow_1,yshadow_1;xshadow_2,yshadow_2;xshadow_3,...
%                yshadow_3;...;xshadow_k,yshadow_k], the kx2 array
%                of the post tip shadow location measurements in the
%                user's local level coordinate system, in meters.
%
%    lpost       The height of the shadow-casting post in meters.
%  
%  Outputs:
%
%    xvecguess   = [delphivecguess;thetasunguess;thetauserguess;...
%                psiuserguess], the first guess of the vector of 
%                unknown paramaters.  delphivecguess(i) is the guessed
%                longitude difference between the sub-solar point
%                and the user point in radians for the measurement
%                time at which the measurements in zhist(i,:) were made
%                for i = 1,...,k, thetasunguess is the guessed latitude
%                of the sub-solar point in radians, thetauserguess is the 
%                guessed latitude of the user in radians, and psiuserguess
%                is the guessed azimuth error of the user's local level
%                coordinate system in radians.  psiuserguess is zero
%                when the user's y axis is aligned with the guessed north.
%                psiuserguess is positive when guessed north is 
%                psiuserguess radians east of (clockwise from) the 
%                user's y axis.  This is a (k+3)x1 vector.
%

%
%  Best fit straight line for the data.
%
   yvsx_linepoly = polyfit(zhist(:,1),zhist(:,2),1);
%
%  Find the nearest point of this curve to the origin and
%  use it to determine the azimuth and the user latitude.  Assume
%  that the sun latitude is zero.
%
   xmin = - (yvsx_linepoly(1,1)*yvsx_linepoly(1,2))/...
                   (1 + yvsx_linepoly(1,1)^2);
   ymin = yvsx_linepoly(1,1)*xmin + yvsx_linepoly(1,2);
   thetasun = 0;
%
%  Use the nearest point of this curve and the assumption that 
%  the sun latitude is zero to figure out the user azimuth error
%  and the user latitude.
%
   if xmin == 0
      psiuser = 0;
   else
      if ymin == 0
         error(['Error in fstgessunlat.m. The initial',...
                 ' azimuth is too far off.'])
      end
      psiuser = atan(xmin/ymin);
   end
   sunelevationatnoon = atan2(lpost,sqrt(xmin^2 + ymin^2));
   thetauser = (0.5*pi - sunelevationatnoon)*sign(ymin);
%
%  Sort through the points to figure out the user longitude differences.
%
   k = size(zhist,1);
   sinthetauser = sin(thetauser);
   costhetauser = cos(thetauser);
   Athetauser = [ sinthetauser, 0, -costhetauser;...
                             0, 1,             0;...
                  costhetauser, 0,  sinthetauser];
   sinpsiuser = sin(psiuser);
   cospsiuser = cos(psiuser);
   Apsiuser = [ -sinpsiuser,  cospsiuser, 0;...
                -cospsiuser, -sinpsiuser, 0;...
                          0,           0, 1];
   Aback = (Apsiuser*Athetauser)';
   xvecguess = [zeros(k,1);thetasun;thetauser;psiuser];
   for j = 1:k
      shatdum = [-zhist(j,:)';lpost];
      shatdum = shatdum*(1/sqrt(shatdum'*shatdum));
      shatdum = Aback*shatdum;
      xvecguess(j,1) = atan2(shatdum(2,1),shatdum(1,1));
   end